// Exercise the actual explorer module with a small DOM adapter and the actual
// shipped IDA Web Worker. This checks the wiring, the rendering and that a
// live run of the compiled model agrees with the committed reference run; it
// checks no browser layout.
import test from "node:test";
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {setTimeout as delay} from "node:timers/promises";
import {Worker} from "node:worker_threads";
import {runInThisContext} from "node:vm";
import {createApp} from "../../docs/pbpk/examples/lobes.mjs";
import {createLobeSolver} from "../../docs/pbpk/examples/lobes-solver.mjs";
import {caseChanges, planSummary} from "../../docs/pbpk/examples/lobes-model.mjs";
import {CONTROLS, inflowChanges} from "../../docs/pbpk/examples/lobes-solver.mjs";

const root = new URL("../../docs/pbpk/examples/", import.meta.url);
const html = await readFile(new URL("lobes.html", root), "utf8");
const data = JSON.parse(await readFile(new URL("lobes-data.json", root), "utf8"));

class Element {
  constructor(tag, attributes = {}) {
    this.tagName = tag;
    this.attributes = attributes;
    this.dataset = Object.fromEntries(Object.entries(attributes).filter(([key]) => key.startsWith("data-")).map(([key, value]) => [key.slice(5), value]));
    this.hidden = "hidden" in attributes;
    this.children = [];
    this.listeners = new Map();
    this.value = attributes.value ?? "";
    this.textContent = "";
    this.className = attributes.class ?? "";
  }
  setAttribute(key, value) {
    this.attributes[key] = String(value);
    if (key === "class") this.className = String(value);
    if (key.startsWith("data-")) this.dataset[key.slice(5)] = String(value);
  }
  getAttribute(key) { return this.attributes[key] ?? null; }
  removeAttribute(key) { delete this.attributes[key]; if (key.startsWith("data-")) delete this.dataset[key.slice(5)]; }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children = children; }
  addEventListener(name, handler) {
    if (!this.listeners.has(name)) this.listeners.set(name, []);
    this.listeners.get(name).push(handler);
  }
  dispatch(name, event = {}) { for (const listener of this.listeners.get(name) ?? []) listener({target: this, preventDefault() {}, ...event}); }
  click() { this.dispatch("click"); }
}

const descendants = element => [element, ...element.children.flatMap(descendants)];

// the named plans sit in a group per operation, so a button of one is found by
// the case of the dataset it carries and not by the row it is in
const presets = byId => byId.get("presets").children.flatMap(group => group.children).filter(node => node.getAttribute("data-case"));
const preset = (byId, name) => presets(byId).find(button => button.getAttribute("data-case") === name);

function stub() {
  const elements = [...html.matchAll(/<([a-z][\w-]*)\s+([^>]*?)>/gi)].map(match => {
    const attributes = Object.fromEntries([...match[2].matchAll(/([\w-]+)(?:="([^"]*)")?/g)].map(a => [a[1], a[2] ?? ""]));
    return new Element(match[1], attributes);
  });
  const byId = new Map(elements.filter(element => element.attributes.id).map(element => [element.attributes.id, element]));
  globalThis.document = {
    getElementById: id => byId.get(id),
    createElement: name => new Element(name),
    createElementNS: (_, name) => new Element(name),
    documentElement: new Element("html"),
  };
  return {
    byId,
    root: {
      querySelector: selector => {
        assert.ok(selector.startsWith("#"), `Unexpected selector: ${selector}`);
        const found = byId.get(selector.slice(1));
        assert.ok(found, `Missing HTML element ${selector}`);
        return found;
      },
    },
  };
}

function fakeSolver(calls) {
  return async () => ({
    configs: {},
    config: () => ({}),
    async simulate(species, controls, changes, options = {}) {
      calls.push({species, controls, changes, options});
      const series = data.species[species].cases.intact.series;
      // a recognizable live result: the remnant falls to half
      return {series: {...series, f_remnant: series.f_remnant.map(() => 0.5)}, provenance: {}};
    },
    close() {},
  });
}

test("the explorer shows a preset, cycles a lobe and runs live when a control moves", async () => {
  const {byId, root: host} = stub();
  const calls = [];
  const app = createApp(host, {data, solverFactory: fakeSolver(calls)});
  await app.render();

  assert.equal(app.state.source, "precomputed");
  assert.equal(byId.get("status").textContent, "Precomputed reference run: the intact liver.");
  const plots = ["plot-masses", "plot-flows", "plot-remnant", "plot-exposure"].map(id => byId.get(id).children[0]);
  for (const plot of plots) assert.equal(plot.tagName, "svg");
  assert.ok(descendants(plots[0]).some(node => node.tagName === "polyline"), "the masses are drawn as curves");
  assert.ok(descendants(plots[1]).some(node => node.tagName === "rect"), "the flows are drawn as bars");

  // the clickable liver has one group per lobe of the species
  const lobes = descendants(byId.get("liver")).filter(node => node.getAttribute("data-lobe"));
  assert.equal(lobes.length, data.species.rat.lobes.length);

  // one click deportalises the lobe, which is a committed case of the dataset
  await app.pick("lll");
  await delay(0);
  assert.equal(app.state.states.lll, "portal");
  assert.deepEqual(caseChanges(app.state.states), {g_po_lll: 0});
  assert.equal(app.state.source, "precomputed", "the dataset carries the single lobe ligation");
  assert.equal(calls.length, 0, "no live run is needed for a committed case at the default controls");

  // a resected lobe leaves the liver and becomes a chip which puts it back
  await app.pick("rl");
  await app.pick("rl");
  await app.pick("rl");
  await delay(0);
  assert.equal(app.state.states.rl, "resected");
  const present = descendants(byId.get("liver")).filter(node => node.getAttribute("data-lobe")).map(node => node.getAttribute("data-lobe"));
  assert.ok(!present.includes("rl"), "the resected lobe holds no area of the liver");
  const chips = byId.get("removed").children;
  assert.deepEqual(chips.map(chip => chip.getAttribute("data-lobe")), ["rl"]);
  chips[0].click();
  await delay(0);
  assert.equal(app.state.states.rl, "intact", "the chip puts the lobe back");
  assert.equal(byId.get("removed").children.length, 0);

  // the presets fill the liver
  const chosen = preset(byId, "hepatectomy70");
  assert.ok(chosen, "the named procedures of the rat are presets");
  chosen.click();
  await delay(0);
  assert.equal(app.state.states.lll, "resected");
  assert.equal(app.state.states.cl, "intact");

  // a slider is a live run of the compiled model
  const before = calls.length;
  const slider = byId.get("sliders").children[0].children[1];
  slider.value = "0.4";
  slider.dispatch("input");
  await delay(0);
  assert.equal(calls.length, before + 1, "the moved control runs the model");
  assert.equal(calls.at(-1).species, "rat");
  assert.equal(calls.at(-1).controls.habrGain, 0.4);
  assert.deepEqual(calls.at(-1).changes, caseChanges(app.state.states));
  assert.equal(app.state.source, "live");
  assert.match(byId.get("status").textContent, /Live run of the compiled model/);

  // the scrubber moves the shown time without running anything
  const ran = calls.length;
  app.setTime(2 * 24 * 60);
  assert.equal(byId.get("clock").textContent, "day 2");
  assert.equal(calls.length, ran);

  // the reset puts the liver and the controls back, so a preset is precomputed again
  byId.get("reset").click();
  await delay(0);
  assert.equal(app.state.states.lll, "intact");
  assert.equal(app.state.controls.habrGain, 0.255);
  assert.equal(app.state.source, "precomputed");

  // every species lays out, including the eight segments of the human
  for (const species of ["mouse", "human"]) {
    await app.setSpecies(species);
    const blocks = descendants(byId.get("liver")).filter(node => node.getAttribute("data-lobe"));
    assert.equal(blocks.length, data.species[species].lobes.length, species);
    assert.equal(app.state.source, "precomputed", species);
  }
  app.close();
});

test("a click during a run is not lost: the page ends on the chosen state", async () => {
  const {byId, root: host} = stub();
  const calls = [];
  let release;
  const gate = new Promise(resolve => { release = resolve; });
  const app = createApp(host, {data, solverFactory: async () => ({
    async simulate(species, controls, changes) {
      calls.push(changes);
      await gate;
      return {series: data.species[species].cases.intact.series, provenance: {}};
    },
    close() {},
  })});
  await app.render();

  // a click while the run of the click before it is still on its way
  const first = app.pick("lml");
  await delay(0);
  assert.equal(calls.length, 1, "the first click runs the model");
  const second = app.pick("rl");
  await delay(0);
  assert.equal(calls.length, 1, "no second run while one is on its way");
  release();
  await Promise.all([first, second]);
  await delay(0);
  assert.equal(calls.length, 2, "the state of the second click is run once the first is back");
  assert.deepEqual(calls.at(-1), {g_po_lml: 0, g_po_rl: 0});
  assert.match(byId.get("status").textContent, /lml PVL, rl PVL/);
  app.close();
});

test("a solver of no species at all is refused", async () => {
  await assert.rejects(createLobeSolver({host: {blobWorkers: () => () => {}}, species: []}), /species/);
});

test("a failing live run keeps the page and says so", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: async () => { throw new Error("no worker here"); }});
  await app.render();
  const slider = byId.get("sliders").children[0].children[1];
  slider.value = "0.5";
  slider.dispatch("input");
  await delay(0);
  assert.match(byId.get("status").textContent, /The live run failed: no worker here/);
  assert.ok(app.state.series, "the last result is kept");
});

test("a live run of the compiled rat agrees with the committed reference run", {timeout: 120000}, async () => {
  stub();
  runInThisContext(await readFile(new URL("solver/sbmlsim-core.js", root), "utf8"));
  const host = globalThis.sbmlsimCore;
  const workers = new Set();
  globalThis.fetch = async url => {
    const bytes = await readFile(new URL(url, root));
    return {
      ok: true,
      status: 200,
      json: async () => JSON.parse(bytes.toString()),
      arrayBuffer: async () => bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength),
    };
  };
  const solver = await createLobeSolver({
    host,
    species: ["rat"],
    createWorker: () => {
      const worker = new Worker(host.WORKER_SOURCE, {eval: true});
      workers.add(worker);
      worker.on("exit", () => workers.delete(worker));
      return worker;
    },
  });
  const changes = data.species.rat.cases.hepatectomy70.changes;
  const result = await solver.simulate("rat", {}, changes);
  const reference = data.species.rat.cases.hepatectomy70.series;
  assert.equal(result.series.time.length, reference.time.length);
  for (const column of ["f_remnant", "Vli_cl", "Qh", "Qpo_cl"]) {
    for (let row = 0; row < reference.time.length; row++) {
      const expected = reference[column][row];
      const actual = result.series[column][row];
      assert.ok(Number.isFinite(actual), `${column}[${row}] is not finite`);
      assert.ok(Math.abs(actual - expected) <= 1e-5 * Math.abs(expected) + 1e-10, `${column}[${row}]: ${actual} != ${expected}`);
    }
  }
  solver.close();
  for (const worker of workers) await worker.terminate();
});

test("the planner selects the lobes of a procedure and the presets fill it", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();

  // the planner offers the intact liver and the three procedures of the library
  const kinds = byId.get("plan-kind").children;
  assert.deepEqual(kinds.map(button => button.getAttribute("data-kind")), ["intact", "portal", "arterial", "resected"]);
  // and one chip per lobe of the species, in the order of the lobe table
  const chips = () => byId.get("plan-lobes").children;
  assert.deepEqual(chips().map(chip => chip.getAttribute("data-lobe")), data.species.rat.lobes.map(lobe => lobe.token));
  for (const chip of chips()) assert.equal(chip.getAttribute("aria-pressed"), "false");

  // clicking a chip puts that lobe into the procedure the planner is set to
  kinds.find(button => button.getAttribute("data-kind") === "resected").click();
  await delay(0);
  chips().find(chip => chip.getAttribute("data-lobe") === "lll").click();
  await delay(0);
  assert.equal(app.state.states.lll, "resected");
  assert.equal(chips().find(chip => chip.getAttribute("data-lobe") === "lll").getAttribute("aria-pressed"), "true");
  // clicking it again takes it back out
  chips().find(chip => chip.getAttribute("data-lobe") === "lll").click();
  await delay(0);
  assert.equal(app.state.states.lll, "intact");

  // a preset sets the procedure and its lobes, and the planner shows both
  preset(byId, "pvl70").click();
  await delay(0);
  assert.equal(app.state.plan.kind, "portal");
  const lit = chips().filter(chip => chip.getAttribute("aria-pressed") === "true").map(chip => chip.getAttribute("data-lobe"));
  assert.deepEqual(lit, ["lll", "lml", "rml"]);
  assert.equal(kinds.find(button => button.getAttribute("data-kind") === "portal").getAttribute("aria-pressed"), "true");

  // the summary of the plan is the mass of the liver it covers
  assert.match(byId.get("plan-summary").textContent, /[Pp]ortal vein ligation of lll, lml, rml/);
  assert.match(byId.get("plan-summary").textContent, /70 %/);
  preset(byId, "hepatectomy90").click();
  await delay(0);
  assert.equal(app.state.plan.kind, "resected");
  assert.match(byId.get("plan-summary").textContent, /10 %/, "the remnant of a 90 % hepatectomy");
  assert.equal(planSummary(data.species.rat.lobes, app.state.states).kind, "resected");
  app.close();
});

test("a setting of the simulation runs the compiled model with what it asks for", async () => {
  const {byId, root: host} = stub();
  const calls = [];
  const app = createApp(host, {data, solverFactory: fakeSolver(calls)});
  await app.render();
  assert.equal(app.state.source, "precomputed");
  assert.equal(calls.length, 0, "the committed run answers the defaults");

  // the page has no run button: a setting of the simulation is a run of its own
  assert.equal(byId.get("run"), undefined, "there is no run button");
  const setting = key => byId.get("sim").children.find(row => row.children[1].getAttribute("data-setting") === key).children[1];
  const days = setting("days");
  days.value = "10";
  days.dispatch("change");
  await delay(0);
  assert.equal(calls.length, 1, "the days which are asked for are run");
  assert.equal(calls.at(-1).options.simulation.days, 10);
  assert.equal(app.state.source, "live");
  // and the time it took is on the page
  assert.match(byId.get("runtime").textContent, /\d+ (ms|s)$/);
  assert.equal(byId.get("spinner").hidden, true, "the spinner is gone once the run is back");

  const dose = setting("dose");
  dose.value = "25";
  dose.dispatch("change");
  await delay(0);
  assert.equal(calls.at(-1).options.simulation.dose, 25);
  assert.equal(calls.at(-1).options.simulation.days, 10, "the settings hold together");

  // the reset puts the settings back, so the committed run answers again
  byId.get("reset").click();
  await delay(0);
  assert.deepEqual(app.state.simulation, {days: 5, points: 121, dose: 10});
  assert.equal(app.state.source, "precomputed");
  app.close();
});

test("a run which is back is played from its first day", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  // the course of a run is played as soon as it is there: a reader sees what
  // the intervention does without asking for it
  assert.equal(app.state.playing, true);
  assert.equal(app.state.time, 0);
  assert.equal(byId.get("play").getAttribute("aria-pressed"), "true");

  app.advance(1000);
  const walked = app.state.time;
  assert.ok(walked > 0, `${walked}`);

  // a plan which is picked plays its own course from the first day again
  preset(byId, "hepatectomy70").click();
  await delay(0);
  assert.equal(app.state.time, 0, "the new course starts at day 0");
  assert.equal(app.state.playing, true);

  // a reader who takes the scrub bar keeps the day they put it on
  const scrub = byId.get("scrub");
  scrub.value = String(2 * 24 * 60);
  scrub.dispatch("input");
  assert.equal(app.state.playing, false, "the scrub bar pauses the course");
  assert.equal(byId.get("clock").textContent, "day 2");
  app.close();
});

test("the spinner is shown while the model runs", async () => {
  const {byId, root: host} = stub();
  let release;
  const gate = new Promise(resolve => { release = resolve; });
  const app = createApp(host, {data, solverFactory: async () => ({
    async simulate(species) { await gate; return {series: data.species[species].cases.intact.series, provenance: {}}; },
    close() {},
  })});
  await app.render();
  assert.equal(byId.get("spinner").hidden, true);
  const running = app.run();
  await delay(0);
  assert.equal(byId.get("spinner").hidden, false, "the spinner says the model is running");
  release();
  await running;
  assert.equal(byId.get("spinner").hidden, true);
  app.close();
});

test("play walks the days of the course and stops at the end", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  const end = data.provenance.end_time;
  // the run which is shown is playing already, so the button pauses it first
  assert.equal(app.state.playing, true);
  byId.get("play").click();
  assert.equal(app.state.playing, false);
  byId.get("play").click();
  assert.equal(app.state.playing, true);
  assert.equal(byId.get("play").getAttribute("aria-pressed"), "true");
  // a frame of the animation moves the clock by the days per second it plays
  app.advance(1000);
  assert.ok(app.state.time > 0 && app.state.time < end, `${app.state.time}`);
  assert.equal(byId.get("scrub").value, String(app.state.time));
  // the clock stops at the last day of the course rather than running past it
  app.advance(1000 * 3600);
  assert.equal(app.state.time, end);
  assert.equal(app.state.playing, false, "it stops at the end");
  assert.equal(byId.get("play").getAttribute("aria-pressed"), "false");
  // playing again from the end starts the course over
  byId.get("play").click();
  assert.equal(app.state.time, 0);
  assert.equal(app.state.playing, true);
  byId.get("play").click();
  assert.equal(app.state.playing, false, "the button pauses as well");
  app.close();
});

test("a lobe is drawn in its own colour and its state is a shade of it", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  const block = token => descendants(byId.get("liver")).find(node => node.getAttribute("data-lobe") === token);
  const part = (token, cls) => descendants(block(token)).find(node => (node.className ?? "").includes(cls));

  // the fill of a block is the colour of its lobe, the same one the curves use
  const colours = Object.fromEntries(data.species.rat.lobes.map((lobe, index) => [lobe.token, `var(--lobe-${index + 1})`]));
  for (const [token, colour] of Object.entries(colours)) {
    assert.match(part(token, "liver-block").getAttribute("style") ?? "", new RegExp(colour.replace(/[()]/g, "\\$&")), token);
  }
  // an intact lobe carries no shading over its colour
  assert.equal(part("lll", "liver-shade").getAttribute("fill"), "none");

  // a ligation shades the box of that lobe and leaves the others alone
  await app.pick("lll");
  await delay(0);
  assert.match(part("lll", "liver-shade").getAttribute("fill"), /^url\(#.*portal.*\)$/);
  assert.equal(part("cl", "liver-shade").getAttribute("fill"), "none");
  // the hepatic artery ligation of the same lobe is a different shading
  await app.pick("lll");
  await delay(0);
  assert.match(part("lll", "liver-shade").getAttribute("fill"), /^url\(#.*arterial.*\)$/);

  // the outline of the whole liver shows what a resection took away
  const frame = () => descendants(byId.get("liver")).find(node => (node.className ?? "").includes("liver-frame"));
  assert.equal(frame().hidden || frame().getAttribute("opacity") === "0", true, "the intact liver needs no outline");
  preset(byId, "hepatectomy70").click();
  await delay(0);
  assert.notEqual(frame().getAttribute("opacity"), "0", "the lost mass is outlined");
  app.close();
});

test("every committed case is recognised, whatever order its changes were written in", async () => {
  // the page matches the plan against the cases of the dataset by its changes;
  // the order the generator wrote them in is not part of what they say
  const {byId, root: host} = stub();
  const calls = [];
  const app = createApp(host, {data, solverFactory: fakeSolver(calls)});
  await app.render();
  for (const species of ["rat", "mouse", "human"]) {
    await app.setSpecies(species);
    for (const name of Object.keys(data.species[species].cases)) {
      preset(byId, name).click();
      await delay(0);
      assert.equal(app.state.source, "precomputed", `${species} ${name} is a committed run`);
    }
  }
  assert.equal(calls.length, 0, "no case of the dataset needs the solver at the default controls");
  app.close();
});

test("a resected lobe keeps a box of its own colour, shaded as resected", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  preset(byId, "hepatectomy70").click();
  await delay(0);
  const chips = byId.get("removed").children;
  assert.deepEqual(chips.map(chip => chip.getAttribute("data-lobe")), ["lll", "lml", "rml"]);
  const box = chip => chip.children.find(child => (child.className ?? "").includes("lobe-dot"));
  for (const chip of chips) {
    const style = box(chip).getAttribute("style") ?? "";
    // the colour of that lobe, the one the diagram and the curves use
    const index = data.species.rat.lobes.findIndex(lobe => lobe.token === chip.getAttribute("data-lobe"));
    assert.ok(style.includes(`var(--lobe-${index + 1})`), `${chip.getAttribute("data-lobe")}: ${style}`);
    // and the shading which says it was taken out
    assert.match(style, /repeating-linear-gradient/);
  }
  // it still puts the lobe back
  chips[0].click();
  await delay(0);
  assert.equal(app.state.states.lll, "intact");
  app.close();
});

test("every lobe hangs between the hilum and the hepatic vein, and a ligation cuts a vessel", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  const paths = () => descendants(byId.get("liver")).filter(node => (node.className ?? "").includes("liver-edge"));
  const of_branch = branch => paths().filter(node => node.className.includes(`liver-edge-${branch}`));

  // three vessels per lobe: the artery and the portal vein in, the hepatic vein out
  for (const branch of ["ha", "po", "hv"]) {
    assert.equal(of_branch(branch).length, data.species.rat.lobes.length, branch);
  }
  // the width of a vessel is its flow, so the hepatic vein of a lobe is the
  // widest of its three (it carries what both inflows brought)
  const width = node => Number(node.getAttribute("stroke-width"));
  assert.ok(width(of_branch("hv")[0]) > width(of_branch("ha")[0]), "the outflow carries both inflows");

  // nothing is cut while the liver is intact
  const cuts = () => descendants(byId.get("liver")).filter(node => node.tagName === "line" && (node.className ?? "").includes("liver-cut"));
  assert.equal(cuts().length, 0);

  // a portal vein ligation cuts that lobe's portal branch, and only it
  await app.pick("lll");
  await delay(0);
  assert.equal(app.state.states.lll, "portal");
  assert.equal(cuts().length, 2, "a cut is two strokes across the vessel");
  assert.ok(cuts().every(node => node.getAttribute("data-branch") === "po_lll"), cuts().map(n => n.getAttribute("data-branch")).join());

  // a hepatic artery ligation cuts the arterial branch instead
  await app.pick("lll");
  await delay(0);
  assert.ok(cuts().every(node => node.getAttribute("data-branch") === "ha_lll"));
  app.close();
});

test("the inflow of a lobe is opened and closed by the controls of the page", () => {
  // the gates of the library are what a ligation closes; the controls scale
  // the same gates, so half the arterial inflow is every open artery at 0.5
  const tokens = ["cl", "lll"];
  const open = {cl: "intact", lll: "intact"};
  assert.deepEqual(inflowChanges(tokens, open, {arterial: 1, portal: 1}), {});
  assert.deepEqual(inflowChanges(tokens, open, {arterial: 0.5, portal: 1}), {g_ha_cl: 0.5, g_ha_lll: 0.5});
  assert.deepEqual(inflowChanges(tokens, open, {arterial: 1, portal: 1.4}), {g_po_cl: 1.4, g_po_lll: 1.4});
  // a ligated branch stays closed whatever the control says
  assert.deepEqual(inflowChanges(tokens, {cl: "portal", lll: "arterial"}, {arterial: 0.5, portal: 1.4}), {g_ha_cl: 0.5, g_po_lll: 1.4});
  // and a resected lobe gets no gate at all
  assert.deepEqual(inflowChanges(tokens, {cl: "resected", lll: "intact"}, {arterial: 0.5, portal: 0.5}), {g_ha_lll: 0.5, g_po_lll: 0.5});
  assert.ok(CONTROLS.arterial && CONTROLS.portal, "the page offers both inflows");
});

test("a moved inflow control runs the model with the gates it asks for", async () => {
  const {byId, root: host} = stub();
  const calls = [];
  const app = createApp(host, {data, solverFactory: fakeSolver(calls)});
  await app.render();
  await app.pick("lll");
  await delay(0);
  const sliders = byId.get("sliders").children;
  const row = sliders.find(entry => entry.children[0].textContent === CONTROLS.arterial.label);
  assert.ok(row, "the arterial inflow is a control of the page");
  row.children[1].value = "0.4";
  row.children[1].dispatch("input");
  await delay(0);
  assert.equal(app.state.controls.arterial, 0.4);
  const changes = calls.at(-1).changes;
  assert.equal(changes.g_po_lll, 0, "the ligation of the lobe stands");
  assert.equal(changes.g_ha_lll, 0.4, "and its artery follows the control");
  assert.equal(changes.g_ha_cl, 0.4);
  app.close();
});

test("the dashboard switches between the volumes and the pressures of a run", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  const titles = () => ["plot-masses", "plot-remnant", "plot-flows", "plot-exposure"]
    .map(id => descendants(byId.get(id)).find(node => (node.className ?? "").includes("cd-title"))?.textContent ?? "");
  const views = byId.get("views").children;
  assert.deepEqual(views.map(button => button.getAttribute("data-view")), ["volume", "pressure"]);
  assert.equal(app.state.view, "volume");
  assert.ok(titles().some(title => /volume of the lobes/.test(title)), titles().join(" | "));

  // the pressure view answers what the pressure network does, on the same run
  views.find(button => button.getAttribute("data-view") === "pressure").click();
  await delay(0);
  assert.equal(app.state.view, "pressure");
  assert.equal(views.find(button => button.getAttribute("data-view") === "pressure").getAttribute("aria-pressed"), "true");
  const shown = titles().join(" | ");
  assert.ok(/pressure/.test(shown), shown);
  assert.ok(/mmHg/.test(shown), shown);
  assert.ok(/perfusion/.test(shown), shown);
  assert.ok(!/volume of the lobes/.test(shown), shown);
  // and the blocks of the liver carry the pressure of the lobe, not its volume
  const value = token => descendants(byId.get("liver"))
    .find(node => node.getAttribute("data-lobe") === token).children
    .find(child => (child.className ?? "").includes("liver-value"))?.textContent ?? "";
  assert.match(value("lll"), /mmHg$/);
  // a block too narrow for a value of nine characters carries none: it would
  // be drawn over its own edge, and the colour and the legend still name it
  assert.equal(value("cl"), "", "the caudate lobe holds 6 % of the mass");

  // no run is needed: it is the same result, read differently
  views.find(button => button.getAttribute("data-view") === "volume").click();
  await delay(0);
  assert.equal(app.state.view, "volume");
  assert.match(value("lll"), /ml$/);
  assert.equal(app.state.source, "precomputed");
  app.close();
});

test("a species is picked by its animal, and the sex is a choice under it", async () => {
  const anatomy = JSON.parse(await readFile(new URL("lobes-anatomy.json", root), "utf8"));
  const {byId, root: host} = stub();
  const app = createApp(host, {data, anatomy, solverFactory: fakeSolver([])});
  await app.render();

  // the three species, each with the animal of the physiology figures on it
  const species = byId.get("species").children;
  assert.deepEqual(species.map(button => button.getAttribute("data-species")), ["rat", "mouse", "human"]);
  assert.equal(species[0].getAttribute("aria-pressed"), "true", "the rat is the species of the page");
  for (const button of species) {
    const icon = button.children.find(child => (child.className ?? "").includes("species-icon"));
    assert.ok(icon, `${button.getAttribute("data-species")} carries its animal`);
    assert.equal(icon.tagName, "svg");
  }

  // the sex is a choice under the species, and offers what that species has
  const sexes = () => byId.get("sex").children;
  assert.deepEqual(sexes().map(button => button.textContent), ["male", "female"]);
  assert.equal(sexes()[1].attributes.disabled, "", "the database has no female rat");
  await app.setSpecies("human");
  assert.equal(sexes()[1].attributes.disabled, undefined, "the database has the female human");
  app.close();
});

test("the named plans are grouped by the operation they are", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();

  // one group per operation the cases of the species cover, in the order of
  // the operations of the planner
  const groups = byId.get("presets").children;
  assert.deepEqual(groups.map(group => group.getAttribute("data-state")), ["intact", "portal", "resected"]);
  assert.deepEqual(groups.map(group => group.children[0].textContent), ["intact", "PVL", "resection"]);

  // a button of a group says what is left of the name of its case: the group
  // above it already says the operation
  const labels = group => group.children.filter(node => node.getAttribute("data-case")).map(node => node.textContent);
  assert.deepEqual(labels(groups[0]), ["whole liver"], "the intact liver is the plan which puts everything back");
  assert.deepEqual(labels(groups[1]), ["70 %", "90 %", "lll"]);
  assert.deepEqual(labels(groups[2]), ["70 %", "90 %", "lll"]);
  for (const button of groups[1].children.filter(node => node.getAttribute("data-case"))) {
    assert.equal(button.getAttribute("data-state"), "portal", "a plan carries the shade of its operation");
  }

  // and the operations of the planner carry the same shade
  const kinds = byId.get("plan-kind").children;
  assert.deepEqual(kinds.map(button => button.getAttribute("data-state")), ["intact", "portal", "arterial", "resected"]);
  app.close();
});

test("the species of the run is drawn, with its liver marked", async () => {
  const {byId, root: host} = stub();
  const anatomy = JSON.parse(await readFile(new URL("lobes-anatomy.json", root), "utf8"));
  const app = createApp(host, {data, anatomy, solverFactory: fakeSolver([])});
  await app.render();
  const figure = () => byId.get("body");
  const svg = () => figure().children.find(child => child.tagName === "svg");
  assert.ok(svg(), "the drawn animal of the species is on the page");
  assert.equal(svg().getAttribute("viewBox"), anatomy.bodies.rat.male.box.join(" "));
  // the liver of the drawing is marked, so a reader sees which organ is split
  const mark = descendants(svg()).find(nodeOf => (nodeOf.className ?? "").includes("body-liver"));
  assert.ok(mark, "the liver of the animal is marked");
  assert.equal(Number(mark.getAttribute("cx")), anatomy.bodies.rat.male.liver[0]);

  // the drawing follows the species
  await app.setSpecies("human");
  assert.equal(svg().getAttribute("viewBox"), anatomy.bodies.human.male.box.join(" "));
  assert.equal(Number(descendants(svg()).find(nodeOf => (nodeOf.className ?? "").includes("body-liver")).getAttribute("cx")), anatomy.bodies.human.male.liver[0]);

  // the licence of the drawings is named on the page
  assert.match(byId.get("provenance").textContent, /CC BY 4\.0/);
  app.close();
});

test("the page works when the drawings are not there", async () => {
  const {byId, root: host} = stub();
  const app = createApp(host, {data, solverFactory: fakeSolver([])});
  await app.render();
  assert.equal(byId.get("body").children.length, 0, "no drawing, no figure");
  assert.equal(byId.get("body").hidden, true);
  assert.ok(byId.get("status").textContent.length > 0, "the rest of the page stands");
  app.close();
});

/** The dataset of the other sex, recognisable by the remnant it reports. */
function femaleData() {
  const female = structuredClone(data);
  female.provenance.sex = "female";
  for (const of_species of Object.values(female.species)) {
    for (const entry of Object.values(of_species.cases)) {
      entry.series.f_remnant = entry.series.f_remnant.map(() => 0.42);
    }
  }
  return female;
}

test("the sex of the individual is a switch, and its dataset is loaded once", async () => {
  const {byId, root: host} = stub();
  const loaded = [];
  const app = createApp(host, {
    data,
    loadData: async sex => { loaded.push(sex); return femaleData(); },
    solverFactory: fakeSolver([]),
  });
  await app.render();
  const buttons = () => byId.get("sex").children;
  assert.deepEqual(buttons().map(child => child.textContent), ["male", "female"]);
  assert.deepEqual(loaded, [], "the sex of the page is loaded with the page");
  // the parameter database has no female rat: the switch says so and refuses
  assert.equal(buttons()[1].getAttribute("disabled"), "");
  assert.match(buttons()[1].getAttribute("title"), /no reference physiology of the female rat/);
  await app.setSex("female");
  assert.equal(app.state.sex, "male", "a sex the species has no physiology of is not taken");
  assert.deepEqual(loaded, []);

  // the human has both, and its female is fetched once
  await app.setSpecies("human");
  assert.ok(!buttons()[1].getAttribute("disabled"), "the female human is offered");
  await app.setSex("female");
  assert.deepEqual(loaded, ["female"]);
  assert.equal(app.state.sex, "female");
  assert.equal(buttons()[1].getAttribute("aria-pressed"), "true");
  // the numbers of the page are the ones of the dataset which was loaded
  assert.ok(
    descendants(byId.get("plot-remnant")).some(node => (node.textContent ?? "").includes("0.42")),
    "the remnant of the female dataset is drawn",
  );

  // both ways, and the dataset of a sex is fetched once
  await app.setSex("male");
  await app.setSex("female");
  assert.deepEqual(loaded, ["female"], "a dataset which is held is not fetched again");

  // a species whose physiology is the male alone takes the page back to the male
  await app.setSpecies("rat");
  assert.equal(app.state.sex, "male");
  assert.match(byId.get("provenance").textContent, /male rat/);
  app.close();
});

test("the drawn body and the solver follow the sex", async () => {
  const {byId, root: host} = stub();
  const anatomy = JSON.parse(await readFile(new URL("lobes-anatomy.json", root), "utf8"));
  const made = [];
  const app = createApp(host, {
    data, anatomy,
    loadData: async () => femaleData(),
    solverFactory: async options => { made.push(options.sex); return (await fakeSolver([])(options)); },
  });
  await app.render();
  await app.setSpecies("human");
  const viewBox = () => byId.get("body").children.find(child => child.tagName === "svg").getAttribute("viewBox");
  assert.equal(viewBox(), anatomy.bodies.human.male.box.join(" "));

  await app.setSex("female");
  assert.equal(viewBox(), anatomy.bodies.human.female.box.join(" "));
  assert.notEqual(anatomy.bodies.human.female.box.join(" "), anatomy.bodies.human.male.box.join(" "));

  // a live run of the female goes to a solver of the models of the female
  await app.run();
  assert.deepEqual(made, ["female"], "the solver is built for the sex which is shown");
  app.close();
});
