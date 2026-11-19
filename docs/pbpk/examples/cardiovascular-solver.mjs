// Browser execution of the compiled SBML heart with sbmlsim-core's shipped
// Web Worker. No cardiovascular equations or alternative integrator live here.

const PARAMETER = {
  hepatic: "CV_fR_hepatic",
  systemic: "CV_fR_systemic",
  portal: "CV_fR_portal",
  renal: "CV_fR_ki",
  heartRate: "CV_HR",
  contractility: "CV_contractility",
};

export function normalizeHeartControls(controls = {}) {
  const values = {hepatic: 1, systemic: 1, portal: 1, renal: 1, heartRate: 75, contractility: 1, ...controls};
  for (const key of Object.keys(PARAMETER)) {
    values[key] = Number(values[key]);
    if (!Number.isFinite(values[key]) || values[key] <= 0) throw new RangeError(`Heart control ${key} must be finite and positive.`);
  }
  if (values.heartRate < 40 || values.heartRate > 150) throw new RangeError("Heart rate must be between 40 and 150/min.");
  return Object.fromEntries(Object.keys(PARAMETER).map(key => [key, values[key]]));
}

function linspace(start, end, steps) {
  return Array.from({length: steps + 1}, (_, i) => start + (end - start) * i / steps);
}

export function heartRequest(controls, config) {
  const values = normalizeHeartControls(controls);
  const beforeStart = config.intervention - config.summary_beats / config.initial_heart_rate;
  // End after an integer number of post-intervention beats. The continuously
  // integrated phase remains continuous because the intervention is at the
  // end of the reference heart's thirtieth beat.
  const finalBeats = Math.floor((config.end - config.intervention) * values.heartRate + 1e-9);
  const end = config.intervention + finalBeats / values.heartRate;
  const afterStart = end - config.summary_beats / values.heartRate;
  const beforeSteps = Math.ceil((config.intervention - beforeStart) * 60 / config.summary_step_seconds);
  const afterSteps = Math.ceil((end - afterStart) * 60 / config.summary_step_seconds);
  const times = [...new Set([
    ...linspace(0, end, 15000),
    ...linspace(beforeStart, config.intervention, beforeSteps),
    ...linspace(afterStart, end, afterSteps),
    ...linspace(afterStart, end, 600),
    beforeStart, config.intervention, afterStart, end,
  ])].sort((a, b) => a - b);
  return {
    contract: "simulation/1",
    start: 0,
    end,
    output: {times},
    preinit: [],
    changes: [{time: config.intervention, values: Object.entries(PARAMETER).map(([key, parameter]) => [parameter, values[key]])}],
    steady_state: null,
    selections: config.columns,
    settings: config.settings,
  };
}

function lowerBound(values, target, column = 0, width = 1) {
  let left = 0, right = values.length / width;
  while (left < right) {
    const middle = Math.floor((left + right) / 2);
    if (values[middle * width + column] < target - 1e-12) left = middle + 1;
    else right = middle;
  }
  return Math.min(left, values.length / width - 1);
}

function summarize(result, start, end, heartRate, before) {
  const width = result.columns.length;
  const time = result.columns.indexOf("time");
  let first = lowerBound(result.values, start, time, width);
  let last = lowerBound(result.values, end, time, width);
  if (!before) {
    while (last + 1 < result.nRows && Math.abs(result.values[(last + 1) * width + time] - end) < 1e-12) last++;
  }
  const index = Object.fromEntries(result.columns.map((name, i) => [name, i]));
  const mean = name => {
    let integral = 0;
    for (let row = first + 1; row <= last; row++) {
      const dt = result.values[row * width + time] - result.values[(row - 1) * width + time];
      integral += 0.5 * dt * (result.values[row * width + index[name]] + result.values[(row - 1) * width + index[name]]);
    }
    return integral / (end - start);
  };
  let systolic = -Infinity, diastolic = Infinity;
  for (let row = first; row <= last; row++) {
    const pressure = result.values[row * width + index.CV_P_ar];
    systolic = Math.max(systolic, pressure); diastolic = Math.min(diastolic, pressure);
  }
  return {
    arterial_systolic: systolic, arterial_diastolic: diastolic,
    mean_arterial_pressure: mean("CV_P_ar"), mean_portal_pressure: mean("CV_P_po"),
    mean_hepatic_pressure: mean("CV_P_li"), mean_cardiac_output: mean("CO"),
    mean_lung_flow: mean("Qlu"), stroke_volume: 1000 * mean("CO") / heartRate,
    heart_rate: heartRate, mean_hepatic_flow: mean("Qh"), mean_portal_flow: mean("Qpo"),
    mean_flows: Object.fromEntries(result.columns.filter(name => name === "CO" || name.startsWith("Q") || name.startsWith("CV_Q_")).map(name => [name, mean(name)])),
  };
}

export function heartCase(result, controls, config) {
  const settings = normalizeHeartControls(controls);
  const width = result.columns.length;
  const timeIndex = result.columns.indexOf("time");
  const end = result.values[(result.nRows - 1) * width + timeIndex];
  const afterStart = end - config.summary_beats / settings.heartRate;
  const table = (times, includeEvent = false) => {
    const rows = new Set(times.map(time => lowerBound(result.values, time, timeIndex, width)));
    if (includeEvent) {
      let event = lowerBound(result.values, config.intervention, timeIndex, width);
      while (event < result.nRows && Math.abs(result.values[event * width + timeIndex] - config.intervention) < 1e-12) rows.add(event++);
    }
    return {columns: result.columns, values: [...rows].sort((a, b) => a - b).map(row => result.columns.map((name, column) => {
      const value = result.values[row * width + column];
      const factor = name === "time" ? 60 : ["CV_V_lv", "CV_V_rv"].includes(name) ? 1000 : 1;
      return Number((value * factor).toFixed(5));
    }))};
  };
  const blood = result.columns.indexOf("CV_V_blood_total");
  const initialBlood = result.values[blood];
  let error = 0;
  for (let row = 0; row < result.nRows; row++) error = Math.max(error, Math.abs(result.values[row * width + blood] - initialBlood));
  return {
    id: "live", label: "Current heart controls", description: "Current resistance, rate and contractility controls, solved with SUNDIALS IDA in a Web Worker.",
    changes: Object.fromEntries(Object.entries(PARAMETER).map(([key, parameter]) => [parameter, settings[key]])),
    controls: settings, heart_rate: settings.heartRate,
    series: table(linspace(0, end, 500), true),
    beats: table(linspace(afterStart, end, 600)),
    summary: {
      before: summarize(result, config.intervention - config.summary_beats / config.initial_heart_rate, config.intervention, config.initial_heart_rate, true),
      after: summarize(result, afterStart, end, settings.heartRate, false),
      maximum_blood_volume_error: error,
    },
    provenance: {engine: "SUNDIALS IDA", host: "sbmlsim-core WebAssembly", sbml_sha256: config.sbml_sha256,
      end_time: end * 60, intervention_time: config.intervention * 60},
  };
}

async function defaultLoadAsset(name) {
  const response = await fetch(new URL(`solver/${name}`, import.meta.url));
  if (!response.ok) throw new Error(`Could not load heart solver ${name}: HTTP ${response.status}.`);
  return new Uint8Array(await response.arrayBuffer());
}

/** A worker of the installed sbmlsim-core host, with its documented protocol. */
export async function createHeartSolver({host = globalThis.sbmlsimCore, loadAsset = defaultLoadAsset, createWorker} = {}) {
  if (!host?.blobWorkers && !createWorker) throw new Error("The sbmlsim-core browser host has not loaded.");
  const [engineBytes, mathlibBytes, modelBytes, descriptionBytes, configBytes] = await Promise.all([
    "engine.wasm", "mathlib.wasm", "model.wasm", "model.json", "config.json",
  ].map(loadAsset));
  const decoder = new TextDecoder();
  const config = JSON.parse(decoder.decode(configBytes));
  const [engine, mathlib, module] = await Promise.all([engineBytes, mathlibBytes, modelBytes].map(bytes => WebAssembly.compile(bytes)));
  const description = decoder.decode(descriptionBytes);
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
    const error = event => stop(new Error(event.message ?? "The IDA simulation worker failed."));
    if (worker.addEventListener) {
      worker.addEventListener("message", receive); worker.addEventListener("error", error);
      detachWorker = () => {
        currentWorker.removeEventListener("message", receive); currentWorker.removeEventListener("error", error);
      };
    } else {
      worker.on("message", receive); worker.on("error", error);
      detachWorker = () => {
        currentWorker.off("message", receive); currentWorker.off("error", error);
      };
    }
    worker.postMessage({type: "init", engine, mathlib});
    worker.postMessage({type: "model", key: "heart", module, description});
  };
  return {
    config,
    async simulate(controls, {signal} = {}) {
      if (closed) throw new Error("The heart solver is closed.");
      signal?.throwIfAborted();
      const normalized = normalizeHeartControls(controls);
      const request = heartRequest(normalized, config);
      if (!worker) spawn();
      const result = await new Promise((resolve, reject) => {
        const id = ++nextId;
        const abort = () => stop(signal.reason ?? new DOMException("Aborted", "AbortError"));
        signal?.addEventListener("abort", abort, {once: true});
        pending.set(id, {resolve, reject, cleanup: () => signal?.removeEventListener("abort", abort)});
        worker.postMessage({type: "simulate", id, key: "heart", request});
      });
      return heartCase(result, normalized, config);
    },
    close() {
      closed = true;
      stop(new DOMException("Closed", "AbortError"));
    },
  };
}
