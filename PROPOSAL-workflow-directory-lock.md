# Proposal: Redis-Based Workflow Directory Lock

## Problem

When a workflow is submitted to the Redis queue, there is no mechanism preventing a second job from being submitted against the same workflow directory while the first is still running. This leads to:

- The second job clearing the `results/` directory while the first is still writing to it
- The `.job_id` file being overwritten, so status polling tracks only the latest job
- Mixed output files from both jobs appearing in the results
- Cross-user file leakage when two sessions share a workspace (via URL query parameter)

## Proposed Solution: Redis-Based Distributed Lock

Use a Redis key as a lightweight lock per workflow directory. The lock is acquired before job submission and released when the job finishes (or fails/times out).

### Why Redis (not filesystem)?

- Already available in the architecture (runs in-container, shared by all workers)
- Atomic operations (`SET NX`) guarantee no race conditions
- Built-in key expiry (`EX`) handles stale locks from crashed workers
- Works across all Streamlit instances and RQ workers without filesystem coordination

## Implementation Plan

### 1. Add lock methods to `QueueManager` (`src/workflow/QueueManager.py`)

Add three methods and one constant:

```python
# Lock TTL should exceed the job timeout to prevent premature expiry,
# but not be so long that a crashed job blocks the directory forever.
LOCK_TTL = 7500  # seconds (job timeout 7200 + 300s buffer)

def acquire_lock(self, workflow_dir: Path) -> bool:
    """
    Attempt to acquire a lock for the given workflow directory.

    Uses Redis SET with NX (only set if not exists) and EX (expiry).
    Returns True if the lock was acquired, False if already held.
    """
    if not self.is_available:
        return True  # No Redis = no lock needed (local mode)

    lock_key = f"lock:workflow:{workflow_dir}"
    # SET NX ensures atomicity: only one caller wins
    return bool(self._redis.set(lock_key, "locked", nx=True, ex=self.LOCK_TTL))

def release_lock(self, workflow_dir: Path) -> None:
    """
    Release the lock for the given workflow directory.
    """
    if not self.is_available:
        return

    lock_key = f"lock:workflow:{workflow_dir}"
    self._redis.delete(lock_key)

def is_locked(self, workflow_dir: Path) -> bool:
    """
    Check if a workflow directory is currently locked.

    Useful for UI feedback (e.g., disabling the start button).
    """
    if not self.is_available:
        return False

    lock_key = f"lock:workflow:{workflow_dir}"
    return bool(self._redis.exists(lock_key))
```

### 2. Acquire lock before job submission (`WorkflowManager._start_workflow_queued`)

In `src/workflow/WorkflowManager.py`, modify `_start_workflow_queued()`:

```python
def _start_workflow_queued(self) -> None:
    """Submit workflow to Redis queue (online mode)"""
    from .tasks import execute_workflow

    # --- NEW: Check lock before submitting ---
    if not self._queue_manager.acquire_lock(self.workflow_dir):
        st.warning("A workflow is already running for this workspace. "
                    "Please wait for it to finish before starting a new one.")
        return
    # --- END NEW ---

    job_id = f"workflow-{self.workflow_dir.name}-{int(time.time())}"

    submitted_id = self._queue_manager.submit_job(
        func=execute_workflow,
        kwargs={
            "workflow_dir": str(self.workflow_dir),
            "workflow_class": self.__class__.__name__,
            "workflow_module": self.__class__.__module__,
        },
        job_id=job_id,
        timeout=7200,
        description=f"Workflow: {self.name}"
    )

    if submitted_id:
        self._queue_manager.store_job_id(self.workflow_dir, submitted_id)
    else:
        # Release lock if submission failed
        self._queue_manager.release_lock(self.workflow_dir)
        st.warning("Queue submission failed, running locally...")
        self._start_workflow_local()
```

### 3. Release lock when the worker task finishes (`tasks.py`)

In `src/workflow/tasks.py`, release the lock in **both** the success and error paths of `execute_workflow()`:

```python
def execute_workflow(workflow_dir: str, workflow_class: str, workflow_module: str) -> dict:
    workflow_path = Path(workflow_dir)

    try:
        # ... existing workflow execution code ...

        return {
            "success": True,
            "workflow_dir": str(workflow_path),
            "message": "Workflow completed successfully"
        }

    except Exception as e:
        # ... existing error handling ...

        return {
            "success": False,
            "workflow_dir": str(workflow_path),
            "error": error_msg
        }

    finally:
        # --- NEW: Always release the lock when done ---
        try:
            from src.workflow.QueueManager import QueueManager
            qm = QueueManager()
            if qm.is_available:
                qm.release_lock(workflow_path)
        except Exception:
            pass  # Lock will expire via TTL if release fails
        # --- END NEW ---
```

The `finally` block guarantees the lock is released regardless of success, failure, or unexpected exceptions. The `LOCK_TTL` expiry acts as a safety net if the worker process crashes entirely (e.g., OOM kill).

### 4. Release lock on job cancellation (`WorkflowManager.stop_workflow`)

In `src/workflow/WorkflowManager.py`, release the lock when a job is canceled:

```python
def stop_workflow(self) -> bool:
    if self._queue_manager and self._queue_manager.is_available:
        job_id = self._queue_manager.load_job_id(self.workflow_dir)
        if job_id:
            success = self._queue_manager.cancel_job(job_id)
            if success:
                self._queue_manager.clear_job_id(self.workflow_dir)
                self._queue_manager.release_lock(self.workflow_dir)  # NEW
                return True

    return self._stop_local_workflow()
```

### 5. UI feedback: disable start button when locked

In `src/workflow/StreamlitUI.py`, check the lock status before rendering the start button. If the directory is locked, show the button as disabled with an informative message:

```python
# In the execution section, before the start button:
is_locked = False
if queue_manager and queue_manager.is_available:
    is_locked = queue_manager.is_locked(workflow_dir)

st.button("Start Workflow", disabled=is_locked)
if is_locked:
    st.info("A workflow is already running. Please wait for it to finish.")
```

## File Summary

| File | Change |
|------|--------|
| `src/workflow/QueueManager.py` | Add `acquire_lock()`, `release_lock()`, `is_locked()`, `LOCK_TTL` |
| `src/workflow/WorkflowManager.py` | Call `acquire_lock()` in `_start_workflow_queued()`, `release_lock()` in `stop_workflow()` |
| `src/workflow/tasks.py` | Call `release_lock()` in `finally` block of `execute_workflow()` |
| `src/workflow/StreamlitUI.py` | Check `is_locked()` to disable start button |

## Edge Cases

| Scenario | Handled by |
|----------|------------|
| Worker crashes mid-execution | `LOCK_TTL` expiry (7500s) auto-releases the lock |
| Job canceled by user | Explicit `release_lock()` in `stop_workflow()` |
| Redis unavailable | `acquire_lock()` returns `True` (falls through, no lock enforced) |
| Queue submission fails after lock acquired | Explicit `release_lock()` in the failure branch |
| Job times out (RQ timeout) | `LOCK_TTL > job_timeout` ensures lock outlives the job; worker `finally` block releases it |

## Limitations / Out of Scope

- **This does not fix local mode.** The lock uses Redis, which is only available in online mode. Local mode uses `multiprocessing.Process` and the PID directory already provides a basic "is running" check. A filesystem-based lock (e.g., `fcntl.flock`) could be added for local mode if needed.
- **This does not prevent URL-based workspace sharing.** Two users accessing the same workspace via URL will still share files. Workspace ownership/authentication is a separate concern. However, the lock does prevent them from running concurrent jobs, which eliminates the most damaging symptom (mixed result files).
