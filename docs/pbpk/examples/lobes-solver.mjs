// Browser execution of the compiled lobed livers with the shipped Web Worker
// of sbmlsim-core. No model of a lobe, no second integrator and no formula of
// the library lives here: this module builds the request of the contract and
// hands the result back as the page's series.

const DOSE = "IVDOSE_drg";
export const CONTROLS = {
  habrGain: {parameter: "habr_gain", label: "buffer response gain", min: 0, max: 0.6, step: 0.005, value: 0.255},
  growth: {parameter: "k_growth", label: "growth constant [1/min]", min: 0, max: 0.002, step: 0.00005, value: 0.0005},
  surgery: {parameter: null, label: "time of the surgery [days]", min: 0, max: 4, step: 0.25, value: 0},
  arterial: {parameter: null, label: "hepatic arterial inflow [x reference]", min: 0, max: 2, step: 0.05, value: 1},
  portal: {parameter: null, label: "portal inflow [x reference]", min: 0, max: 2, step: 0.05, value: 1},
};

/** The gates of the lobes an inflow control asks for.
 *
 * A branch of a lobe is opened and closed by the gate `pkpdsim.pbpk.lobes`
 * writes (`g_ha_<lobe>`, `g_po_<lobe>`), which a ligation sets to zero; the
 * controls scale the same gate of every branch which is still open, so the
 * page changes the inflow of a liver without a second model of a vessel. A
 * branch the plan ligated keeps its zero and is left to `caseChanges`, and a
 * resected lobe has no open branch at all.
 */
export function inflowChanges(tokens, states = {}, controls = {}) {
  const values = normalizeControls(controls);
  const changes = {};
  for (const token of tokens) {
    const state = states[token] ?? "intact";
    if (state === "resected") continue;
    if (state !== "arterial" && values.arterial !== CONTROLS.arterial.value) changes[`g_ha_${token}`] = values.arterial;
    if (state !== "portal" && values.portal !== CONTROLS.portal.value) changes[`g_po_${token}`] = values.portal;
  }
  return changes;
}

/** The controls of the page with every value checked, the defaults filled in. */
export function normalizeControls(controls = {}) {
  const values = {};
  for (const [key, control] of Object.entries(CONTROLS)) {
    const value = Number(controls[key] ?? control.value);
    if (!Number.isFinite(value)) throw new RangeError(`The control ${key} must be a finite number.`);
    if (value < control.min || value > control.max) throw new RangeError(`The control ${key} must be between ${control.min} and ${control.max}.`);
    values[key] = value;
  }
  return values;
}

function linspace(start, end, steps) {
  return Array.from({length: steps + 1}, (_, i) => start + ((end - start) * i) / steps);
}

/** The request of one run: the dose and the controls before the initialization,
 * the intervention at its time, and the times of the output.
 *
 * A surgery at time 0 is applied before the initialization, which is how a
 * study which simulates only the days after an operation sets it up; a later
 * one is a change of the running model, which is a surgery at a time point.
 */
export function lobeRequest(controls, changes, config) {
  const values = normalizeControls(controls);
  const surgery = values.surgery * 24 * 60;
  const settings = [
    [CONTROLS.habrGain.parameter, values.habrGain],
    [CONTROLS.growth.parameter, values.growth],
  ];
  const intervention = Object.entries(changes ?? {});
  const times = [...new Set([
    ...linspace(config.start, config.end, Math.max(config.points - 1, 1)),
    ...(surgery > config.start && surgery < config.end ? [surgery, surgery * (1 + 1e-9)] : []),
  ])].sort((a, b) => a - b);
  return {
    contract: "simulation/1",
    start: config.start,
    end: config.end,
    output: {times},
    preinit: [[DOSE, config.dose], ...settings, ...(surgery <= config.start ? intervention : [])],
    changes: surgery > config.start && intervention.length > 0 ? [{time: surgery, values: intervention}] : [],
    steady_state: null,
    selections: config.columns,
    settings: config.settings,
  };
}

/** The result of the worker as the series of the page, one array per column. */
export function lobeSeries(result) {
  const width = result.columns.length;
  const series = Object.fromEntries(result.columns.map(name => [name, new Array(result.nRows)]));
  for (let row = 0; row < result.nRows; row++) {
    for (let column = 0; column < width; column++) series[result.columns[column]][row] = result.values[row * width + column];
  }
  return series;
}

export const RUNTIME_ASSETS = ["engine.wasm", "mathlib.wasm", "sbmlsim-core.js"];

/** Where an asset of a species lives: the model beside it, the runtime shared.
 *
 * The engine and the host script are the ones of the cardiovascular examples,
 * `solver/`, so the bytes of the engine are served once for both pages
 * (scripts/lobes_examples.py writes only the model).
 */
export function assetPath(species, name, sex = "male") {
  if (RUNTIME_ASSETS.includes(name)) return `solver/${name}`;
  const folder = sex === "male" ? "" : `${sex}/`;
  return `${folder}solver_lobes/${species}/${name}`;
}

async function defaultLoadAsset(species, name, sex) {
  const response = await fetch(new URL(assetPath(species, name, sex), import.meta.url));
  if (!response.ok) throw new Error(`Could not load ${name} of the ${species}: HTTP ${response.status}.`);
  return new Uint8Array(await response.arrayBuffer());
}

/** A worker of the installed sbmlsim-core host which runs the lobed liver of a species.
 *
 * One worker holds the compiled models of every species which was asked for,
 * under the key of the species, so switching species does not recompile. The
 * models of a sex are its own, so a worker is a worker of one sex.
 */
export async function createLobeSolver({host = globalThis.sbmlsimCore, loadAsset = defaultLoadAsset, createWorker, species = [], sex = "male"} = {}) {
  if (species.length === 0) throw new RangeError("A solver of the lobes needs the species whose models it runs.");
  if (!host?.blobWorkers && !createWorker) throw new Error("The sbmlsim-core browser host has not loaded.");
  const decoder = new TextDecoder();
  const configs = {};
  const models = {};
  for (const name of species) {
    const configBytes = await loadAsset(name, "config.json", sex);
    configs[name] = JSON.parse(decoder.decode(configBytes));
    const [modelBytes, descriptionBytes] = await Promise.all([
      loadAsset(name, "model.wasm", sex),
      loadAsset(name, "model.json", sex),
    ]);
    models[name] = {module: await WebAssembly.compile(modelBytes), description: decoder.decode(descriptionBytes)};
  }
  const [engineBytes, mathlibBytes] = await Promise.all([
    loadAsset(species[0], "engine.wasm"),
    loadAsset(species[0], "mathlib.wasm"),
  ]);
  const [engine, mathlib] = await Promise.all([engineBytes, mathlibBytes].map(bytes => WebAssembly.compile(bytes)));
  const factory = createWorker ?? host.blobWorkers();
  let worker, detachWorker, nextId = 0, closed = false;
  const pending = new Map();
  const stop = error => {
    detachWorker?.(); detachWorker = null;
    worker?.terminate(); worker = null;
    for (const job of pending.values()) { job.cleanup(); job.reject(error); }
    pending.clear();
  };
  const spawn = () => {
    worker = factory();
    const currentWorker = worker;
    const receive = event => {
      const message = event?.data ?? event;
      const job = pending.get(message.id);
      if (!job) return;
      pending.delete(message.id); job.cleanup();
      if (message.type === "error") job.reject(Object.assign(new Error(message.error.message), message.error));
      else job.resolve(message);
    };
    const error = event => stop(new Error(event.message ?? "The IDA simulation of the lobes failed."));
    if (worker.addEventListener) {
      worker.addEventListener("message", receive); worker.addEventListener("error", error);
      detachWorker = () => { currentWorker.removeEventListener("message", receive); currentWorker.removeEventListener("error", error); };
    } else {
      worker.on("message", receive); worker.on("error", error);
      detachWorker = () => { currentWorker.off("message", receive); currentWorker.off("error", error); };
    }
    worker.postMessage({type: "init", engine, mathlib});
    for (const [name, model] of Object.entries(models)) worker.postMessage({type: "model", key: name, module: model.module, description: model.description});
  };
  return {
    configs,
    config(name) { return configs[name]; },
    async simulate(name, controls, changes, {signal} = {}) {
      if (closed) throw new Error("The solver of the lobes is closed.");
      if (!configs[name]) throw new RangeError(`The ${name} has no compiled model on this page.`);
      signal?.throwIfAborted();
      const request = lobeRequest(controls, changes, configs[name]);
      if (!worker) spawn();
      const result = await new Promise((resolve, reject) => {
        const id = ++nextId;
        const abort = () => stop(signal.reason ?? new DOMException("Aborted", "AbortError"));
        signal?.addEventListener("abort", abort, {once: true});
        pending.set(id, {resolve, reject, cleanup: () => signal?.removeEventListener("abort", abort)});
        worker.postMessage({type: "simulate", id, key: name, request});
      });
      return {
        series: lobeSeries(result),
        provenance: {engine: "SUNDIALS IDA", host: "sbmlsim-core WebAssembly", sbml_sha256: configs[name].sbml_sha256, species: name},
      };
    },
    close() { closed = true; stop(new DOMException("Closed", "AbortError")); },
  };
}
