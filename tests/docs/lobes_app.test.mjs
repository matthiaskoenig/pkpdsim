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
import {caseChanges} from "../../docs/pbpk/examples/lobes-model.mjs";

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
  setAttribute(key, value) { this.attributes[key] = String(value); if (key.startsWith("data-")) this.dataset[key.slice(5)] = String(value); }
  getAttribute(key) { return this.attributes[key] ?? null; }
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
    async simulate(species, controls, changes) {
      calls.push({species, controls, changes});
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
  const preset = byId.get("presets").children.find(button => button.textContent === "hepatectomy70");
  assert.ok(preset, "the named procedures of the rat are presets");
  preset.click();
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
