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
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
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
  if (files.length === 0) throw new Error("Export: no files written");
};
