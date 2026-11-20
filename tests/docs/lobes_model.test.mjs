import test from "node:test";
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {LOBE_SLOTS, STATES, caseChanges, csv, duration, exposure, lobeColors, planSummary, sampleAt, stateOfCase, trim} from "../../docs/pbpk/examples/lobes-model.mjs";
import {blockLayout, cutMark, flowWidth} from "../../docs/pbpk/examples/lobes-schematic.mjs";

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

test("a lobe keeps its colour whatever else is shown", () => {
  // colour follows the lobe, not its rank in what is plotted: a resected lobe
  // leaving the diagram must not repaint the ones which stay
  const colours = lobeColors(rat.lobes);
  assert.deepEqual(Object.keys(colours), tokens);
  assert.deepEqual([...new Set(Object.values(colours))].length, tokens.length);
  for (const [index, token] of tokens.entries()) {
    assert.equal(colours[token], `var(--lobe-${index + 1})`);
  }
  // the eight segments of the human each get their own slot of the palette
  const human = lobeColors(data.species.human.lobes);
  assert.equal(Object.keys(human).length, 8);
  assert.equal([...new Set(Object.values(human))].length, 8);
  // past the slots of the palette a lobe is not given a cycled hue
  const many = lobeColors([...data.species.human.lobes, {token: "x9", name: "ninth"}]);
  assert.equal(many.x9, "var(--lobe-other)");
  assert.equal(LOBE_SLOTS, 8);
});

test("the summary of a plan is the mass it takes and the mass which is left", () => {
  // the share of the liver mass a procedure covers, which is what the names of
  // the procedures of the literature count (pvl70 is 70 % of the mass)
  const intact = Object.fromEntries(tokens.map(token => [token, "intact"]));
  const none = planSummary(rat.lobes, intact);
  close(none.affected, 0);
  close(none.remnant, 1);
  assert.equal(none.kind, undefined);
  const pvl70 = planSummary(rat.lobes, {...intact, lll: "portal", lml: "portal", rml: "portal"});
  close(pvl70.affected, 0.7);
  close(pvl70.remnant, 1, 1e-9);
  assert.equal(pvl70.kind, "portal");
  assert.deepEqual(pvl70.lobes, ["lll", "lml", "rml"]);
  // a resection takes the mass out of the liver, a ligation leaves it in place
  const h70 = planSummary(rat.lobes, {...intact, lll: "resected", lml: "resected", rml: "resected"});
  close(h70.affected, 0.7);
  close(h70.remnant, 0.3);
  assert.equal(h70.kind, "resected");
  // two kinds at once is a mixed plan and says so
  const mixed = planSummary(rat.lobes, {...intact, lll: "portal", rl: "resected"});
  assert.equal(mixed.kind, "mixed");
  close(mixed.remnant, 0.8);
});

test("the named procedures of the dataset are the shares their names say", () => {
  for (const [name, entry] of Object.entries(rat.cases)) {
    const summary = planSummary(rat.lobes, stateOfCase(entry, tokens));
    const digits = name.match(/(\d+)$/);
    if (!digits) continue;
    const share = Number(digits[1]) / 100;
    close(summary.affected, share, 1e-9);
    if (name.startsWith("hepatectomy")) close(summary.remnant, 1 - share, 1e-9);
  }
});

test("a run time is shown in the unit a reader reads", () => {
  assert.equal(duration(0), "0 ms");
  assert.equal(duration(7.4), "7 ms");
  assert.equal(duration(942), "942 ms");
  assert.equal(duration(1500), "1.5 s");
  assert.equal(duration(62000), "62 s");
  assert.equal(duration(Number.NaN), "-");
});

test("the area of the liver is the mass which is there, not the frame", () => {
  // the total area of the blocks is the size of the liver: after a resection
  // of 70 % of the mass the blocks hold 30 % of the area they held intact
  const frame = {width: 400, height: 200};
  const shares = Object.fromEntries(rat.lobes.map(lobe => [lobe.token, lobe.mass_fraction]));
  const full = blockLayout(rat.lobes, {...frame, shares, scale: 1});
  const area = blocks => blocks.reduce((sum, block) => sum + block.width * block.height, 0);
  close(area(full), frame.width * frame.height, 1e-6);
  const third = blockLayout(rat.lobes, {...frame, shares, scale: 0.3});
  close(area(third) / area(full), 0.3, 1e-9);
  // every block stays inside the frame and the small liver is centred in it
  for (const block of third) {
    assert.ok(block.x >= -1e-9 && block.y >= -1e-9, "no block leaves the frame");
    assert.ok(block.x + block.width <= frame.width + 1e-9);
    assert.ok(block.y + block.height <= frame.height + 1e-9);
  }
  const left = Math.min(...third.map(block => block.x));
  const right = Math.max(...third.map(block => block.x + block.width));
  close((left + (frame.width - right)) / 2, left, 1e-6);
  // and it sits on the hilum, so the vessels keep their length
  close(Math.max(...third.map(block => block.y + block.height)), frame.height, 1e-6);
  // the shares within the liver are unchanged by the scale
  for (const block of third) {
    const was = full.find(entry => entry.token === block.token);
    close(block.width * block.height / area(third), was.width * was.height / area(full), 1e-9);
  }
  // without a scale the blocks fill the frame, as they did before
  close(area(blockLayout(rat.lobes, {...frame, shares})), frame.width * frame.height, 1e-6);
});

test("a lobe keeps its place in the liver and only grows or shrinks", () => {
  // the drawing is read over time: a lobe which atrophies while its neighbour
  // grows must stay where it was, so the rows and the order of the blocks come
  // from the lobe table and never from the shares of the moment
  const frame = {width: 400, height: 240};
  const reference = blockLayout(rat.lobes, frame);
  const order = reference.map(block => block.token);
  const place = blocks => Object.fromEntries(blocks.map(block => [block.token, {x: block.x, y: block.y}]));
  const was = place(reference);

  // the smallest lobe grows past the largest one, which is what a remnant does
  const skewed = {cl: 0.62, lll: 0.1, lml: 0.1, rl: 0.1, rml: 0.08};
  const grown = blockLayout(rat.lobes, {...frame, shares: skewed});
  assert.deepEqual(grown.map(block => block.token), order, "the blocks keep their order");
  const now = place(grown);
  // the lobes stand in columns: a lobe keeps its column, so the left to right
  // order never swaps, and keeps the lobes above and below it in that column
  const column = (blocks, token) => blocks.find(block => block.token === token).x;
  for (const a of order) {
    for (const b of order) {
      if (a === b) continue;
      const together = column(reference, a) === column(reference, b);
      assert.equal(together, column(grown, a) === column(grown, b), `${a} and ${b} changed column`);
      if (together) assert.equal(Math.sign(now[a].y - now[b].y), Math.sign(was[a].y - was[b].y), `${a} and ${b} swapped places in their column`);
      else assert.equal(Math.sign(now[a].x - now[b].x), Math.sign(was[a].x - was[b].x), `${a} and ${b} swapped sides`);
    }
  }
  // and the areas are the shares of the moment, filling the frame
  const area = blocks => blocks.reduce((sum, block) => sum + block.width * block.height, 0);
  close(area(grown) / (frame.width * frame.height), 1, 1e-6);
  for (const block of grown) close(block.width * block.height / area(grown), skewed[block.token], 1e-9);

  // a resected lobe leaves and the ones which stay keep their order
  const gone = blockLayout(rat.lobes, {...frame, shares: {...skewed, lml: 0}});
  assert.deepEqual(gone.map(block => block.token), order.filter(token => token !== "lml"));
  close(area(gone) / (frame.width * frame.height), 1, 1e-6);
});

test("a ligated vessel is cut across its middle", () => {
  // a straight edge from (0,0) to (0,100): the cut is two strokes across it,
  // at its middle and square to it
  const cut = cutMark([0, 0], [0, 33], [0, 67], [0, 100]);
  assert.equal(cut.length, 2);
  for (const [a, b] of cut) {
    close((a[1] + b[1]) / 2, 50, 4);
    close(a[1], b[1], 1e-9);
  }
  // and it is put where the vessel is free of its neighbours, by a fraction
  // of the curve: at three quarters it is three quarters of the way along
  for (const [a, b] of cutMark([0, 0], [0, 33], [0, 67], [0, 100], {at: 0.75})) {
    close((a[1] + b[1]) / 2, 75, 4);
  }
  for (const [a, b] of cut) {
    close(a[1], b[1], 1e-9, "the stroke is square to the vessel");
    assert.ok(Math.abs(a[0] - b[0]) > 6, "the stroke is wide enough to read");
  }
  assert.ok(Math.abs(cut[0][0][1] - cut[1][0][1]) > 2, "the two strokes are apart");
  // a horizontal vessel is cut vertically
  const across = cutMark([0, 0], [33, 0], [67, 0], [100, 0]);
  for (const [a, b] of across) close(a[0], b[0], 1e-9);
  // a vessel of no length gets no cut rather than a NaN
  assert.deepEqual(cutMark([5, 5], [5, 5], [5, 5], [5, 5]), []);
});
