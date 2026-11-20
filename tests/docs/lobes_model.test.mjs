import test from "node:test";
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {STATES, caseChanges, csv, exposure, sampleAt, stateOfCase, trim} from "../../docs/pbpk/examples/lobes-model.mjs";
import {blockLayout, flowWidth} from "../../docs/pbpk/examples/lobes-schematic.mjs";

const data = JSON.parse(await readFile(new URL("../../docs/pbpk/examples/lobes-data.json", import.meta.url), "utf8"));
const rat = data.species.rat;
const tokens = rat.lobes.map(lobe => lobe.token);
const close = (actual, expected, tolerance = 1e-9) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);

test("the four states of a lobe are the interventions of the library", () => {
  assert.deepEqual(STATES.map(state => state.id), ["intact", "portal", "arterial", "resected"]);
  for (const state of STATES) assert.ok(state.label.length > 0 && state.description.length > 0);
});

test("the changes of a state are the parameters the library changes", () => {
  assert.deepEqual(caseChanges({lll: "intact", cl: "intact"}), {});
  assert.deepEqual(caseChanges({lll: "portal"}), {g_po_lll: 0});
  assert.deepEqual(caseChanges({lll: "arterial"}), {g_ha_lll: 0});
  assert.deepEqual(caseChanges({lll: "resected"}), {g_ha_lll: 0, g_po_lll: 0, f_present_lll: 0});
  assert.deepEqual(caseChanges({lll: "portal", rl: "resected"}), {g_po_lll: 0, g_ha_rl: 0, g_po_rl: 0, f_present_rl: 0});
});

test("the state of every committed case is read back from its changes", () => {
  for (const [name, entry] of Object.entries(rat.cases)) {
    const state = stateOfCase(entry, tokens);
    assert.deepEqual(Object.keys(state).sort(), [...tokens].sort(), name);
    for (const token of entry.lobes) {
      assert.equal(state[token], entry.kind === "Resection" ? "resected" : "portal", `${name} ${token}`);
    }
    for (const token of tokens.filter(t => !entry.lobes.includes(t))) assert.equal(state[token], "intact", `${name} ${token}`);
    // what the UI builds from the state is what the generator applied
    assert.deepEqual(caseChanges(state), entry.changes, name);
  }
});

test("a case of the dataset carries the course of every lobe", () => {
  const series = rat.cases.hepatectomy70.series;
  assert.equal(series.Vli_cl.length, series.time.length);
  assert.ok(Number(series.Vli_cl.at(-1)) > 2 * Number(series.Vli_cl[0]), "the remnant of a 70 % hepatectomy grows");
  assert.deepEqual(sampleAt(undefined, 0), {});
});

test("the exposure of a lobe is the trapezoid area of its outflow", () => {
  close(exposure([0, 1, 2], [0, 2, 0]), 2);
  close(exposure([0, 10], [1, 1]), 10);
  close(exposure([5], [3]), 0);
  const area = exposure(rat.cases.intact.series.time, rat.cases.intact.series.Cout_li_cl_drg);
  assert.ok(area > 0, "the intact liver clears the dose");
});

test("a sample is the row at or before a time", () => {
  const time = rat.cases.intact.series.time;
  const sample = sampleAt(rat.cases.intact.series, 0);
  assert.equal(sample.time, time[0]);
  assert.equal(sampleAt(rat.cases.intact.series, 1e9).time, time.at(-1));
  assert.ok(sampleAt(rat.cases.intact.series, 1440).time <= 1440);
  assert.ok(Number.isFinite(sampleAt(rat.cases.intact.series, 1440).Qh));
});

test("the csv of a case carries its columns and its rows", () => {
  const text = csv(rat.cases.pvl70.series, ["time", "Vli_lll", "Qpo_lll"]);
  const lines = text.trim().split("\n");
  assert.equal(lines[0], "time,Vli_lll,Qpo_lll");
  assert.equal(lines.length, rat.cases.pvl70.series.time.length + 1);
  assert.equal(lines[1].split(",").length, 3);
  assert.throws(() => csv(rat.cases.pvl70.series, ["nope"]), /nope/);
});

test("a number is trimmed for a label and never shown as a false zero", () => {
  assert.equal(trim(0), "0");
  assert.equal(trim(1.23456), "1.235");
  assert.equal(trim(0.000123), "1.23e-4");
  assert.equal(trim(Number.NaN), "-");
  assert.equal(trim(Infinity), "-");
});

test("the blocks of a liver have the area of the shares of its lobes", () => {
  const layout = blockLayout(rat.lobes, {width: 300, height: 200});
  assert.equal(layout.length, rat.lobes.length);
  const total = layout.reduce((sum, block) => sum + block.width * block.height, 0);
  close(total / (300 * 200), 1, 1e-6);
  for (const block of layout) {
    assert.ok(block.width > 0 && block.height > 0, block.token);
    assert.ok(block.x >= -1e-9 && block.y >= -1e-9 && block.x + block.width <= 300 + 1e-9 && block.y + block.height <= 200 + 1e-9, block.token);
  }
  // the largest lobe gets the largest block
  const largest = [...layout].sort((a, b) => b.width * b.height - a.width * a.height)[0];
  assert.equal(largest.token, rat.lobes.reduce((a, b) => (a.mass_fraction > b.mass_fraction ? a : b)).token);
});

test("a lobe with no mass left is not a block of the liver", () => {
  const shares = Object.fromEntries(rat.lobes.map(lobe => [lobe.token, lobe.mass_fraction]));
  shares.lll = 0;
  const layout = blockLayout(rat.lobes, {width: 300, height: 200, shares});
  assert.deepEqual(layout.map(block => block.token).sort(), tokens.filter(t => t !== "lll").sort());
  const total = layout.reduce((sum, block) => sum + block.width * block.height, 0);
  close(total / (300 * 200), 1, 1e-6);
});

test("a liver whose every lobe is gone has no block and no stray geometry", () => {
  const shares = Object.fromEntries(rat.lobes.map(lobe => [lobe.token, 0]));
  assert.deepEqual(blockLayout(rat.lobes, {width: 300, height: 200, shares}), []);
});

test("the eight segments of the human and the five lobes of the mouse lay out too", () => {
  for (const species of ["human", "mouse"]) {
    const layout = blockLayout(data.species[species].lobes, {width: 320, height: 240});
    assert.equal(layout.length, data.species[species].lobes.length);
    const total = layout.reduce((sum, block) => sum + block.width * block.height, 0);
    close(total / (320 * 240), 1, 1e-6);
  }
});

test("an edge grows with the square root of its flow and never vanishes", () => {
  close(flowWidth(0), 0.7);
  assert.ok(flowWidth(1) > flowWidth(0.25));
  close(flowWidth(4) - flowWidth(1), flowWidth(1) - flowWidth(0), 1e-9);
  close(flowWidth(-1), flowWidth(1));
  close(flowWidth(Number.NaN), 0.7);
});
