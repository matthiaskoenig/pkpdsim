import test from "node:test";
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {column, csv, defaultControls, lastBeat, sampleAt, solvePressure} from "../../docs/pbpk/examples/cardiovascular-model.mjs";
import {flowWidth, schematicValues} from "../../docs/pbpk/examples/cardiovascular-schematic.mjs";

const data = JSON.parse(await readFile(new URL("../../docs/pbpk/examples/cardiovascular-data.json", import.meta.url), "utf8"));
const reference = data.pressure;
const close = (actual, expected, tolerance = 1e-10) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);

test("steady browser calculation recovers the exported SBML reference", () => {
  const result = solvePressure(reference, defaultControls(reference));
  close(result.ppo, reference.portal_pressure);
  close(result.pli, reference.hepatic_pressure);
  close(result.co, reference.cardiac_output);
  for (const organ of result.flows) close(organ.flow, reference.organs.find(item => item.id === organ.id).flow);
});

test("slider range conserves volume at both junctions", () => {
  for (const hepatic of [1, 2, 4, 8]) for (const systemic of [0.5, 1, 2]) {
    for (const portal of [1, 6]) for (const renal of [0.5, 3]) {
      const controls = {...defaultControls(reference), hepatic, systemic, portal, renal};
      const result = solvePressure(reference, controls);
      close(result.portalResidual, 0);
      close(result.hepaticResidual, 0);
      close(result.co, result.systemicFlow + result.qh);
      assert.ok(result.flows.every(organ => organ.flow > 0));
      assert.ok(controls.arterial > result.ppo && result.ppo > result.pli && result.pli > controls.venous);
    }
  }
});

test("hepatic resistance raises portal pressure and reduces portal flow", () => {
  const baseline = solvePressure(reference, defaultControls(reference));
  const changed = solvePressure(reference, {...defaultControls(reference), hepatic: 4});
  assert.ok(changed.ppo > baseline.ppo);
  assert.ok(changed.pli > baseline.pli);
  assert.ok(changed.qpo < baseline.qpo);
  close(changed.flows.find(organ => organ.id === "ki").flow, baseline.flows.find(organ => organ.id === "ki").flow);
});

test("renal resistance changes its own flow without changing other branches", () => {
  const baseline = solvePressure(reference, defaultControls(reference));
  const changed = solvePressure(reference, {...defaultControls(reference), renal: 2});
  for (const organ of changed.flows) {
    const expected = baseline.flows.find(item => item.id === organ.id).flow / (organ.id === "ki" ? 2 : 1);
    close(organ.flow, expected);
  }
});

test("invalid boundaries and resistances are rejected", () => {
  for (const change of [{hepatic: 0}, {systemic: -1}, {arterial: 5, venous: 10}, {portal: NaN}]) {
    assert.throws(() => solvePressure(reference, {...defaultControls(reference), ...change}), RangeError);
  }
});

test("heartbeat views select complete beats and preserve dimensional CSV headers", () => {
  for (const example of data.cases) {
    const selected = lastBeat(example.beats, example.heart_rate);
    const times = column(selected, "time");
    assert.ok(times.at(-1) - times[0] <= 60 / example.heart_rate + 1e-8);
    assert.ok(times.at(-1) - times[0] > 0.99 * 60 / example.heart_rate);
    assert.ok(column(selected, "CV_V_lv").every(value => value > 0));
    assert.match(csv(selected, data.provenance.units), /"CV_V_lv \[ml\]"/);
    assert.match(csv(selected, data.provenance.units), /"time \[s\]"/);
  }
  assert.throws(() => column(data.cases[0].series, "missing"), /Missing reference column/);
});

test("replay samples clamp at endpoints and respect the post-intervention event row", () => {
  const trace = {columns: ["Q", "time"], values: [[1, 0], [1, 1], [3, 1], [3, 2]]};
  assert.deepEqual(sampleAt(trace, -1), {Q: 1, time: 0});
  assert.deepEqual(sampleAt(trace, 0.5), {Q: 1, time: 0.5});
  assert.deepEqual(sampleAt(trace, 1), {Q: 3, time: 1});
  assert.deepEqual(sampleAt(trace, 1.5), {Q: 3, time: 1.5});
  assert.deepEqual(sampleAt(trace, 3), {Q: 3, time: 2});
  const ramp = {columns: ["time", "CV_P_ar"], values: [[0, 80], [1, 120]]};
  close(sampleAt(ramp, 0.25).CV_P_ar, 90);
});

test("schematic widths use one monotonic magnitude scale and separate portal storage flows", () => {
  const widths = [0, 1, 5, 20, 100].map(flow => flowWidth(flow));
  assert.ok(widths.every(Number.isFinite));
  assert.ok(widths.every((width, index) => !index || width > widths[index - 1]));
  close(flowWidth(5, 2) - flowWidth(0), 2 * (flowWidth(5) - flowWidth(0)));
  close(flowWidth(-5), flowWidth(5));
  close(flowWidth(NaN), flowWidth(0));
  close(flowWidth(4) - flowWidth(0), 2 * (flowWidth(1) - flowWidth(0)));
  const values = schematicValues({Qhe: 0.3, Qbr: 0.8, Qmu: 1.2, Qki: 1.3,
    Qst: 0.1, Qsi: 0.7, Qla: 0.3, Qsp: 0.2, Qpa: 0.1, Qpo: 1.1, Qha: 0.4});
  close(values.systemic, 2.3);
  close(values.splanchnic, 1.4);
  close(values.Qpo, 1.1);
});
