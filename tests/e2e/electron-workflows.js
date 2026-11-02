// UmetaFlow scenario for electron/scripts/e2e.js: loads the bundled example
// data set (four mzML files), runs the UmetaFlow workflow with its default
// parameters inside the Electron window and checks that the results page
// shows the feature matrix and draws its plot, then exports the result files
// into a folder.
const fs = require("fs");
const os = require("os");
const path = require("path");

const RUN_TIMEOUT = 60 * 60 * 1000;

module.exports = async (t) => {
  await t.open("file_upload");
  await t.clickTab("Example Data");
  await t.clickButton("Load Example Data");
  await t.assertNoException("File Upload");

  await t.open("umetaflow_run");
  await t.clickButton("Start Workflow");
  try {
    await t.waitForText("Workflow completed successfully.", {
      timeout: RUN_TIMEOUT,
      failOn: ["Errors occurred, check log file."],
    });
  } catch (err) {
    printWorkflowLog(t);
    throw err;
  }
  await t.shot("workflow-finished");

  await t.open("umetaflow_results");
  await t.waitForText("Feature Matrix", { timeout: 3 * 60 * 1000 });
  const started = Date.now();
  let plots = 0;
  while (Date.now() - started < 3 * 60 * 1000) {
    await t.assertNoException("Results");
    plots = await t.page.locator(".js-plotly-plot").count();
    if (plots > 0) break;
    await t.page.waitForTimeout(2000);
  }
  t.log(`Results: ${plots} plots`);
  await t.shot("results");
  if (plots === 0) throw new Error("Results: no plots drawn");

  const exportDir = fs.mkdtempSync(path.join(os.tmpdir(), "umetaflow-export-"));
  await t.page.getByRole("button", { name: "📤 Export" }).first().click();
  const folder = t.page.getByLabel("Export folder");
  await folder.fill(exportDir);
  await folder.press("Enter");
  await t.settle();
  await t.page.getByRole("button", { name: "📤 Export" }).last().click();
  await t.waitForText("Exported", { timeout: 3 * 60 * 1000 });
  const exported = path.join(exportDir, "UmetaFlow-results");
  const files = fs.existsSync(exported) ? fs.readdirSync(exported) : [];
  t.log(`Export: ${files.length} entries in ${exported}`);
  await t.shot("export");
  if (!files.includes("meta-value-template.tsv")) throw new Error("Export: meta-value-template.tsv missing");
  if (!files.some((f) => f.endsWith(".parquet"))) throw new Error("Export: feature matrix parquet missing");
};

// The workflow's own log is in the workspace, not in the Streamlit log the CI
// job prints, so a failed run shows it here.
function printWorkflowLog(t) {
  const appName = "UmetaFlow";
  const roots = [
    path.join(os.homedir(), ".config", appName),
    path.join(os.homedir(), "Library", "Application Support", appName),
    path.join(process.env.APPDATA || path.join(os.homedir(), "AppData", "Roaming"), appName),
  ];
  for (const root of roots) {
    const log = path.join(root, "workspaces", "default", "umetaflow", "logs", "all.log");
    if (!fs.existsSync(log)) continue;
    const lines = fs.readFileSync(log, "utf8").split(/\r?\n/);
    t.log(`Last lines of ${log}:`);
    console.log(lines.slice(-150).join("\n"));
    return;
  }
  t.log("No workflow log found");
}
