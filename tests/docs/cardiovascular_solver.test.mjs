// The actual shipped sbmlsim-core Web Worker and WASM runtime execute in a
// Node worker thread here. This is numerical integration, not a DOM mock.
import test from "node:test";
import assert from "node:assert/strict";
import {readFile, writeFile} from "node:fs/promises";
import {Worker} from "node:worker_threads";
import {runInThisContext} from "node:vm";
import {createHeartSolver, heartRequest, normalizeHeartControls} from "../../docs/pbpk/examples/cardiovascular-solver.mjs";

const directory = new URL("../../docs/pbpk/examples/solver/", import.meta.url);
runInThisContext(await readFile(new URL("sbmlsim-core.js", directory), "utf8"));
const sex = process.env.PKPDSIM_HEART_TEST_SEX ?? "male";
const modelDirectory = sex === "female" ? new URL("../female/solver/", directory) : directory;
const loadAsset = name => readFile(new URL(name, ["engine.wasm", "mathlib.wasm"].includes(name) ? directory : modelDirectory));
const config = JSON.parse(await readFile(new URL("config.json", modelDirectory), "utf8"));
const host = globalThis.sbmlsimCore;

test("heart requests preserve shared baseline, complete beats and requested controls", () => {
  const controls = normalizeHeartControls({hepatic: 4, systemic: 0.8, portal: 1.7, renal: 2, heartRate: 87, contractility: 0.9});
  const request = heartRequest(controls, config);
  assert.equal(request.changes[0].time, 0.4);
  assert.equal(request.preinit.length, 0);
  assert.equal(request.changes[0].values.find(([name]) => name === "CV_fR_ki")[1], 2);
  assert.ok(Math.abs((request.end - 0.4) * 87 - Math.round((request.end - 0.4) * 87)) < 1e-10);
  assert.ok(request.output.times.every((time, index, times) => !index || time >= times[index - 1]));
  assert.throws(() => normalizeHeartControls({hepatic: 0}), RangeError);
  assert.throws(() => normalizeHeartControls({heartRate: 1}), RangeError);
});

test("the browser worker solves arbitrary simultaneous resistance changes with IDA", async () => {
  const solver = await createHeartSolver({host, loadAsset, createWorker: () => new Worker(host.WORKER_SOURCE, {eval: true})});
  try {
    const baseline = await solver.simulate({});
    const controls = {hepatic: 4, systemic: 0.8, portal: 1.7, renal: 2, heartRate: 87, contractility: 0.9};
    const changed = await solver.simulate(controls);
    assert.equal(changed.provenance.engine, "SUNDIALS IDA");
    assert.equal(changed.heart_rate, 87);
    assert.equal(changed.beats.values.length, 601);
    assert.ok(changed.summary.after.mean_portal_pressure > baseline.summary.after.mean_portal_pressure);
    assert.ok(changed.summary.after.mean_cardiac_output > baseline.summary.after.mean_cardiac_output);
    assert.ok(changed.summary.maximum_blood_volume_error < 1e-9);
    assert.ok(Math.abs(changed.summary.after.mean_cardiac_output - changed.summary.after.mean_lung_flow) < 1e-3);
    assert.ok(changed.series.columns.includes("CV_Q_tv"));
    assert.ok(changed.series.columns.includes("CV_V_rv"));
    assert.ok(changed.summary.after.mean_flows.Qki > 0);
    const columns = Object.fromEntries(changed.beats.columns.map((name, i) => [name, i]));
    for (const row of changed.beats.values) {
      assert.ok(row.every(Number.isFinite));
      assert.ok(Math.abs(row[columns.Qha] + row[columns.Qpo] - row[columns.Qh]) < 3e-5);
    }
    if (process.env.PKPDSIM_HEART_TEST_OUTPUT) {
      await writeFile(process.env.PKPDSIM_HEART_TEST_OUTPUT, JSON.stringify({sex, controls, baseline, changed}));
    }
    const controller = new AbortController();
    const pending = solver.simulate({hepatic: 8}, {signal: controller.signal});
    controller.abort();
    await assert.rejects(pending, {name: "AbortError"});
    const recovered = await solver.simulate({});
    assert.ok(Math.abs(recovered.summary.after.mean_cardiac_output - baseline.summary.after.mean_cardiac_output) < 1e-8);
  } finally { solver.close(); }
  await assert.rejects(solver.simulate({}), /closed/);
});
