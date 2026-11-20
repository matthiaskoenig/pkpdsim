// Exercise the actual app module with a small DOM adapter and the actual
// shipped IDA Web Worker. This checks wiring and rendering, not browser layout.
import test from "node:test";
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {setImmediate as nextTurn, setTimeout as delay} from "node:timers/promises";
import {Worker} from "node:worker_threads";
import {runInThisContext} from "node:vm";
import {flowWidth} from "../../docs/pbpk/examples/cardiovascular-schematic.mjs";
import {solvePressure} from "../../docs/pbpk/examples/cardiovascular-model.mjs";

const root = new URL("../../docs/pbpk/examples/", import.meta.url);
const html = await readFile(new URL("cardiovascular.html", root), "utf8");
const data = JSON.parse(await readFile(new URL("cardiovascular-data.json", root), "utf8"));

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
    this.clientWidth = 640;
    this.height = 1600;
  }
  setAttribute(key, value) { this.attributes[key] = String(value); }
  getAttribute(key) { return this.attributes[key] ?? null; }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children = children; }
  addEventListener(name, handler) {
    if (!this.listeners.has(name)) this.listeners.set(name, []);
    this.listeners.get(name).push(handler);
  }
  dispatch(name, fields = {}) { for (const listener of this.listeners.get(name) ?? []) listener({target: this, ...fields}); }
  click() { this.dispatch("click"); }
  getBoundingClientRect() { return {height: this.height, width: this.clientWidth, left: 0, top: 0}; }
}

function descendants(element) {
  return [element, ...element.children.flatMap(descendants)];
}

async function waitUntil(predicate, message, timeout = 15000) {
  const deadline = Date.now() + timeout;
  while (!predicate()) {
    assert.ok(Date.now() < deadline, message);
    await delay(20);
  }
}

test("the app shares sliders, solves combined controls, replays flow widths and recovers from loading errors", {timeout: 40000}, async () => {
  const elements = [...html.matchAll(/<([a-z][\w-]*)\s+([^>]*?)>/gi)].map(match => {
    const attributes = Object.fromEntries([...match[2].matchAll(/([\w-]+)(?:="([^"]*)")?/g)].map(attribute => [attribute[1], attribute[2] ?? ""]));
    return new Element(match[1], attributes);
  });
  const byId = new Map(elements.filter(element => element.attributes.id).map(element => [element.attributes.id, element]));
  const select = selector => {
    if (selector.startsWith("[data-")) {
      const name = selector.slice(1, -1);
      return elements.filter(element => name in element.attributes);
    }
    throw new Error(`Unexpected selector: ${selector}`);
  };
  globalThis.document = {
    getElementById: id => { assert.ok(byId.has(id), `Missing HTML element ${id}`); return byId.get(id); },
    createElementNS: (_, name) => new Element(name),
    createElement: name => new Element(name),
    querySelectorAll: select,
    documentElement: new Element("html"),
  };
  byId.get("flow-speed").value = "0.5";
  const messages = [];
  const windowEvents = new Element("window");
  globalThis.window = {location: {origin: "http://localhost"}, parent: {
    document: {body: new Element("body", {"data-md-color-scheme": "default"})},
    postMessage: value => messages.push(value),
  }, addEventListener: (...args) => windowEvents.addEventListener(...args)};
  globalThis.getComputedStyle = () => ({getPropertyValue: () => "#007c78"});
  globalThis.matchMedia = query => ({matches: query.includes("reduced-motion"), addEventListener() {}});
  globalThis.ResizeObserver = class { observe() {} };
  globalThis.MutationObserver = class { observe() {} };
  const frames = new Map();
  let nextFrame = 0;
  globalThis.requestAnimationFrame = callback => { frames.set(++nextFrame, callback); return nextFrame; };
  globalThis.cancelAnimationFrame = id => frames.delete(id);
  const advanceFrame = time => {
    const ready = [...frames.values()];
    frames.clear();
    for (const callback of ready) callback(time);
  };

  runInThisContext(await readFile(new URL("solver/sbmlsim-core.js", root), "utf8"));
  const host = globalThis.sbmlsimCore;
  const workers = new Set();
  let workerCount = 0;
  globalThis.sbmlsimCore = {...host, blobWorkers: () => () => {
    const worker = new Worker(host.WORKER_SOURCE, {eval: true});
    workers.add(worker); workerCount++;
    worker.on("exit", () => workers.delete(worker));
    return worker;
  }};
  let shouldFail = true;
  const assets = [];
  globalThis.fetch = async url => {
    if (shouldFail) return {ok: false, status: 503};
    const asset = new URL(url, root);
    assets.push(asset.pathname);
    const bytes = await readFile(asset);
    return {ok: true, json: async () => JSON.parse(bytes.toString()),
      arrayBuffer: async () => bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)};
  };
  const errors = [];
  const originalError = console.error;
  console.error = (...args) => errors.push(args);
  const mode = name => select("[data-mode]").find(button => button.dataset.mode === name).click();
  const edges = () => descendants(byId.get("circulation-schematic")).filter(element => element.getAttribute("data-edge"));
  const edge = name => {
    const group = edges().find(element => element.getAttribute("data-edge") === name);
    assert.ok(group, `Missing flow edge ${name}`);
    return group.children.find(child => child.getAttribute("data-flow") !== null);
  };
  const scrub = time => { byId.get("flow-time").value = String(time); byId.get("flow-time").dispatch("input"); };
  const assertFiniteCharts = () => {
    for (const chart of ["arterial-chart", "portal-chart", "pv-chart", "flow-chart"]) {
      const svg = byId.get(chart).children[0];
      assert.equal(svg.tagName, "svg");
      const paths = svg.children.filter(child => child.tagName === "path");
      assert.equal(paths.length, 2);
      for (const path of paths) assert.doesNotMatch(path.attributes.d, /NaN|Infinity/);
      for (const child of svg.children.filter(item => item.tagName === "line" || item.tagName === "circle")) {
        for (const coordinate of ["x1", "x2", "y1", "y2", "cx", "cy"]) {
          if (child.getAttribute(coordinate) !== null) assert.ok(Number.isFinite(Number(child.getAttribute(coordinate))));
        }
      }
    }
  };
  try {
    await import(new URL("cardiovascular.mjs", root));
    await nextTurn();
    assert.equal(byId.get("load-error").hidden, false);
    assert.equal(byId.get("app").getAttribute("aria-busy"), "false");
    assert.equal(errors.length, 1);
    shouldFail = false;
    byId.get("retry").click();
    await waitUntil(() => !byId.get("workspace").hidden, "Reference data did not load");
    assert.equal(byId.get("load-error").hidden, true);
    assert.equal(byId.get("metric-0").textContent, "10.0");
    assert.equal(byId.get("heart-case").children.length, 9);
    assert.equal(edges().length, 34);

    select("[data-preset]").find(button => button.dataset.preset === "hepatic").click();
    assert.ok(Number(byId.get("metric-0").textContent) > 15);
    byId.get("reset").click();
    assert.equal(byId.get("metric-0").textContent, "10.0");
    byId.get("hepatic").value = "8";
    byId.get("hepatic").dispatch("input");
    assert.ok(Number(byId.get("metric-0").textContent) > 25);
    assert.equal(byId.get("organ-flows").children.length, data.pressure.organs.length + 1);

    mode("heart");
    assert.equal(byId.get("pressure-controls").hidden, false);
    assert.equal(byId.get("pressure-boundaries").hidden, true);
    assert.equal(byId.get("pressure-results").hidden, true);
    assert.equal(byId.get("heart-results").hidden, false);
    assert.equal(Number(byId.get("hepatic").value), 8);
    assert.equal(edges().length, 38);
    assert.equal(frames.size, 0, "Reduced motion must prevent automatic replay");
    const parameters = {hepatic: "CV_fR_hepatic", systemic: "CV_fR_systemic", portal: "CV_fR_portal", renal: "CV_fR_ki", "heart-rate": "CV_HR", contractility: "CV_contractility"};
    for (const example of data.cases) {
      byId.get("heart-case").value = example.id;
      byId.get("heart-case").dispatch("change");
      assert.equal(byId.get("case-description").textContent, example.description);
      for (const [id, parameter] of Object.entries(parameters)) assert.equal(Number(byId.get(id).value), example.changes[parameter] ?? (id === "heart-rate" ? 75 : 1));
      for (const interval of select("[data-window]")) { interval.click(); assertFiniteCharts(); }
    }
    assert.equal(byId.get("metric-3").textContent, "100");
    assert.equal(byId.get("metric-label-0").textContent, "Mean portal pressure");

    // Inspect a closed and an open aortic valve in the same solved beat.
    byId.get("heart-case").value = "baseline";
    byId.get("heart-case").dispatch("change");
    const baseline = data.cases.find(example => example.id === "baseline");
    const timeColumn = baseline.beats.columns.indexOf("time");
    const coColumn = baseline.beats.columns.indexOf("CO");
    const rows = baseline.beats.values;
    const peak = rows.reduce((current, row) => row[coColumn] > current[coColumn] ? row : current);
    const closed = rows.find(row => row[coColumn] === 0);
    assert.ok(closed, "The reference trajectory should include a closed aortic valve");
    scrub(closed[timeColumn] - rows[0][timeColumn]);
    const closedWidth = Number(edge("Aortic valve").getAttribute("stroke-width"));
    assert.equal(Number(edge("Aortic valve").getAttribute("data-flow")), 0);
    assert.equal(closedWidth, flowWidth(0));
    scrub(peak[timeColumn] - rows[0][timeColumn]);
    const valve = edge("Aortic valve");
    assert.ok(Math.abs(Number(valve.getAttribute("data-flow")) - peak[coColumn]) < 1e-4);
    assert.ok(Number(valve.getAttribute("stroke-width")) > closedWidth);
    for (const group of edges()) {
      const line = group.children.find(child => child.getAttribute("data-flow") !== null);
      const value = line.getAttribute("data-flow");
      assert.notEqual(value, "", "Every schematic edge needs a simulated flow");
      const flow = Number(value);
      assert.ok(Math.abs(Number(line.getAttribute("stroke-width")) - flowWidth(flow)) < 1e-12);
    }
    assert.equal(edge("Renal arterial supply").getAttribute("stroke-width"), edge("Renal venous return").getAttribute("stroke-width"));
    assertFiniteCharts();
    scrub(0);
    byId.get("play-flow").click();
    assert.equal(byId.get("play-flow").getAttribute("aria-pressed"), "true");
    assert.equal(frames.size, 1);
    advanceFrame(1000); advanceFrame(1200);
    const slowTime = Number(byId.get("flow-time").value);
    assert.ok(Math.abs(slowTime - 0.1) < 1e-10);
    byId.get("flow-speed").value = "2";
    byId.get("flow-speed").dispatch("change");
    advanceFrame(1400); advanceFrame(1500);
    assert.ok(Math.abs(Number(byId.get("flow-time").value) - slowTime - 0.2) < 1e-10);
    assert.equal(workerCount, 0, "Replay and playback speed must not start a new simulation");
    assertFiniteCharts();
    byId.get("play-flow").click();
    assert.equal(frames.size, 0);
    assert.equal(byId.get("play-flow").getAttribute("aria-pressed"), "false");

    // Run a genuine worker for a combination absent from the preset cache.
    const combined = {hepatic: 4, systemic: 0.8, portal: 1.7, renal: 2, "heart-rate": 87, contractility: 0.9};
    for (const [id, value] of Object.entries(combined)) {
      byId.get(id).value = String(value);
      byId.get(id).dispatch("input");
    }
    assert.equal(byId.get("heart-case").value, "custom");
    await waitUntil(() => byId.get("solver-status").textContent.startsWith("Current settings solved"), `Live IDA solve failed: ${byId.get("solver-status").textContent}`);
    assert.equal(workerCount, 1);
    assert.ok(assets.some(asset => asset.endsWith("/solver/engine.wasm")));
    assert.equal(byId.get("metric-3").textContent, "87");
    assert.ok(Number(byId.get("metric-0").textContent) > baseline.summary.after.mean_portal_pressure);
    assert.ok(Number(byId.get("metric-1").textContent) > baseline.summary.after.mean_cardiac_output);
    assert.ok(Number(byId.get("conservation").textContent.match(/error: ([\deE.+-]+)/)[1]) < 1e-9);
    const flow = name => Number(edge(name).getAttribute("data-flow"));
    assert.ok(Math.abs(flow("Hepatic arterial supply") + flow("Portal vein to sinusoids") - flow("Hepatic venous return")) < 5e-5);
    assertFiniteCharts();

    let downloaded;
    const originalCreate = URL.createObjectURL;
    URL.createObjectURL = blob => { downloaded = blob; return originalCreate(blob); };
    try {
      byId.get("download").click();
      assert.match(await downloaded.text(), /CV_V_lv \[ml\]/);
    } finally { URL.createObjectURL = originalCreate; }
    byId.get("app").height = 1100;
    mode("pressure");
    assert.equal(edges().length, 34);
    assert.equal(byId.get("heart-results").hidden, true);
    assert.equal(byId.get("pressure-boundaries").hidden, false);
    assert.equal(byId.get("metric-label-0").textContent, "Portal pressure");
    for (const [id, value] of Object.entries(combined)) assert.equal(Number(byId.get(id).value), value);
    const expected = solvePressure(data.pressure, {...combined, arterial: 95, venous: 5});
    assert.ok(Math.abs(Number(byId.get("metric-0").textContent) - expected.ppo) < 0.051);
    assert.equal(messages.at(-1).height, 1100);
    mode("heart");
    assert.equal(workerCount, 1, "Switching back should reuse the solved slider combination");
    // Display scaling changes every vessel and legend without solving the model.
    scrub(0.15);
    const hit = descendants(byId.get("arterial-chart")).find(item => item.getAttribute("class") === "inspection-hit");
    hit.dispatch("pointerdown", {clientX: 400, pointerId: 1});
    assert.ok(Number(byId.get("flow-time").value) > 0.5, "Clicking the chart should scrub the displayed instant");
    assert.equal(frames.size, 0);
    const previousFlow = flow("Aortic valve");
    byId.get("flow-scale").value = "2";
    byId.get("flow-scale").dispatch("input");
    assert.equal(workerCount, 1);
    assert.equal(flow("Aortic valve"), previousFlow);
    assert.equal(Number(edge("Aortic valve").getAttribute("stroke-width")), flowWidth(previousFlow, 2));
    for (const chart of ["arterial-chart", "portal-chart", "flow-chart", "pv-chart"]) {
      const marker = descendants(byId.get(chart)).find(item => item.getAttribute("data-cursor-point") === chart);
      assert.ok(marker, `Missing synchronized point in ${chart}`);
      assert.equal(marker.getAttribute("r"), "5.5");
    }
    const maleBody = descendants(byId.get("circulation-schematic")).find(item => item.getAttribute("class") === "cv-body").innerHTML;
    byId.get("sex").value = "female";
    byId.get("sex").dispatch("change");
    await waitUntil(() => !byId.get("workspace").hidden, "Female physiology did not load");
    assert.equal(byId.get("individual-label").textContent, "Female reference");
    assert.equal(Number(byId.get("hepatic").value), combined.hepatic);
    await waitUntil(() => byId.get("solver-status").textContent.startsWith("Current settings solved"), "Female IDA solve failed");
    assert.equal(workerCount, 2);
    assert.ok(assets.some(asset => asset.endsWith("/female/solver/model.wasm")));
    assert.notEqual(descendants(byId.get("circulation-schematic")).find(item => item.getAttribute("class") === "cv-body").innerHTML, maleBody);
    const femaleOutput = byId.get("metric-1").textContent;
    byId.get("sex").value = "male";
    byId.get("sex").dispatch("change");
    await waitUntil(() => !byId.get("workspace").hidden, "Male physiology did not reload");
    assert.equal(workerCount, 2, "A cached male result must not be replaced by the female result");
    assert.notEqual(byId.get("metric-1").textContent, femaleOutput);
    assertFiniteCharts();
    window.innerWidth = 1186; window.parent.innerHeight = 900;
    byId.get("flow-scale").dispatch("input");
    assert.equal(messages.at(-1).height, 780, "A wide documentation frame should fill its available viewport height");
    window.innerWidth = 390;
    byId.get("flow-scale").dispatch("input");
    assert.equal(messages.at(-1).height, 1100, "Narrow frames should follow the natural content height");
    assert.equal(errors.length, 1, "Only the intentional first loading failure should be logged");
  } finally {
    windowEvents.dispatch("pagehide");
    await waitUntil(() => workers.size === 0, "Page teardown did not close the IDA worker", 3000);
    assert.equal(frames.size, 0);
    console.error = originalError;
  }
});
