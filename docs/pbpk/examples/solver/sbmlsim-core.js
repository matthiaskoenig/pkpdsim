// sbmlsim-core 0.24.0, the inline build of its JavaScript host (js/inline.mjs, js/build.mjs).
// Generated, do not edit. MIT AND BSD-3-Clause AND LGPL-2.1-or-later, see LICENSE and THIRD_PARTY.md.
var sbmlsimCore = (() => {
  var __defProp = Object.defineProperty;
  var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
  var __getOwnPropNames = Object.getOwnPropertyNames;
  var __hasOwnProp = Object.prototype.hasOwnProperty;
  var __export = (target, all) => {
    for (var name in all)
      __defProp(target, name, { get: all[name], enumerable: true });
  };
  var __copyProps = (to, from, except, desc) => {
    if (from && typeof from === "object" || typeof from === "function") {
      for (let key of __getOwnPropNames(from))
        if (!__hasOwnProp.call(to, key) && key !== except)
          __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
    }
    return to;
  };
  var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

  // inline.mjs
  var inline_exports = {};
  __export(inline_exports, {
    BUNDLE_CONTRACT: () => BUNDLE_CONTRACT,
    CORE_VERSION: () => CORE_VERSION,
    ENGINE: () => ENGINE,
    WORKER_SOURCE: () => WORKER_SOURCE,
    blobWorkers: () => blobWorkers,
    createInlineEndpoint: () => createInlineEndpoint
  });

  // loader.mjs
  var MODEL_CONTRACT = "model/2";
  var MAX_MEMORY_PAGES = 32768;
  var ENGINE_SECTION = "sbmlsim.engine";
  var ENGINE_OUTPUT_FORMAT = 2;
  var WASI = "wasi_snapshot_preview1";
  var ENGINE_EXPORTS = [
    "engine_alloc",
    "engine_free",
    "engine_load",
    "engine_simulate",
    "engine_output_ptr",
    "engine_output_len",
    "engine_stats_ptr"
  ];
  var STATS_NAMES = [
    "steps",
    "rhs_evals",
    "jacobian_evals",
    "jacobian_rhs_evals",
    "jacobian_fallbacks",
    "linear_setups",
    "error_test_failures",
    "nonlinear_iterations",
    "nonlinear_convergence_failures",
    "root_evals",
    "integrations",
    "probe_rhs_evals",
    "steady_state_rhs_evals",
    "linear_iterations",
    "linear_convergence_failures",
    "jacobian_vector_evals",
    "preconditioner_evals",
    "preconditioner_solves",
    "sensitivity_rhs_evals",
    "sensitivity_model_rhs_evals",
    "sensitivity_jac_p_evals"
  ];
  var LITTLE_ENDIAN = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
  var EngineError = class extends Error {
    constructor(kind, message, { flag, time } = {}) {
      super(`${kind}: ${message}`);
      this.name = "EngineError";
      this.kind = kind;
      this.detail = message;
      if (flag !== void 0) this.flag = flag;
      if (time !== void 0) this.time = time === null ? NaN : time;
    }
  };
  var UnsupportedError = class extends EngineError {
    constructor(message) {
      super("unsupported", message);
      this.name = "UnsupportedError";
    }
  };
  var WasiExit = class extends Error {
    constructor(code) {
      super(`engine.wasm called proc_exit(${code})`);
      this.name = "WasiExit";
      this.code = code;
    }
  };
  async function toModule(source, what) {
    if (source instanceof WebAssembly.Module) return source;
    try {
      return await WebAssembly.compile(source);
    } catch (error) {
      throw new Error(`${what} does not compile: ${error.message}`);
    }
  }
  function memoryMinimum(engineModule) {
    const sections = WebAssembly.Module.customSections(engineModule, ENGINE_SECTION);
    if (sections.length !== 1) {
      throw new Error(
        `engine.wasm has no custom section ${ENGINE_SECTION} with the minimum of its memory and the version of its output (the build script of sbmlsim-host adds it)`
      );
    }
    const { memory_minimum_pages: pages, output_format: format } = JSON.parse(
      new TextDecoder().decode(sections[0])
    );
    if (format === void 0) {
      throw new Error(
        `engine.wasm declares no output format in its section ${ENGINE_SECTION}, this loader reads ${ENGINE_OUTPUT_FORMAT}`
      );
    }
    if (format !== ENGINE_OUTPUT_FORMAT) {
      throw new Error(
        `engine.wasm writes the output format ${format}, this loader reads ${ENGINE_OUTPUT_FORMAT}`
      );
    }
    if (!Number.isInteger(pages) || pages < 0 || pages > MAX_MEMORY_PAGES) {
      throw new Error(`engine.wasm declares a memory of ${pages} pages`);
    }
    return pages;
  }
  function wasiImports(engineModule, getMemory, out) {
    let calls = 0;
    const decoders = /* @__PURE__ */ new Map();
    const lines = /* @__PURE__ */ new Map();
    const view = () => new DataView(getMemory().buffer);
    const write = (fd, bytes) => {
      if (!decoders.has(fd)) decoders.set(fd, new TextDecoder());
      let text = (lines.get(fd) ?? "") + decoders.get(fd).decode(bytes, { stream: true });
      let newline;
      while ((newline = text.indexOf("\n")) >= 0) {
        (fd === 2 ? out.error : out.log)(text.slice(0, newline));
        text = text.slice(newline + 1);
      }
      lines.set(fd, text);
    };
    const provided = {
      fd_write(fd, iovs, iovsLen, nwritten) {
        calls += 1;
        if (fd !== 1 && fd !== 2) return 8;
        const memory = view();
        let total = 0;
        for (let i = 0; i < iovsLen; i += 1) {
          const ptr = memory.getUint32(iovs + 8 * i, true);
          const len = memory.getUint32(iovs + 8 * i + 4, true);
          write(fd, new Uint8Array(getMemory().buffer, ptr, len));
          total += len;
        }
        memory.setUint32(nwritten, total, true);
        return 0;
      },
      proc_exit(code) {
        calls += 1;
        throw new WasiExit(code);
      },
      environ_get() {
        calls += 1;
        return 0;
      },
      environ_sizes_get(count, size) {
        calls += 1;
        const memory = view();
        memory.setUint32(count, 0, true);
        memory.setUint32(size, 0, true);
        return 0;
      },
      clock_time_get() {
        calls += 1;
        throw new Error("engine.wasm called WASI clock_time_get: the engine reads no clock");
      },
      random_get() {
        calls += 1;
        throw new Error("engine.wasm called WASI random_get: the engine draws no random bytes");
      }
    };
    const imports = {};
    for (const { module, name, kind } of WebAssembly.Module.imports(engineModule)) {
      if (module !== WASI || kind !== "function") continue;
      imports[name] = provided[name] ?? (() => {
        calls += 1;
        throw new Error(`engine.wasm called WASI ${name}, which the loader does not provide`);
      });
    }
    return { imports, calls: () => calls };
  }
  function checkRuntime(engineModule, mathlibModule) {
    let memory = false;
    for (const { module, name, kind } of WebAssembly.Module.imports(engineModule)) {
      if (module === "env" && name === "memory" && kind === "memory") memory = true;
      else if ((module === "model" || module === WASI) && kind === "function") continue;
      else throw new Error(`engine.wasm imports ${module}.${name}`);
    }
    if (!memory) throw new Error("engine.wasm imports no memory");
    const exports = new Set(WebAssembly.Module.exports(engineModule).map((e) => e.name));
    for (const name of ENGINE_EXPORTS) {
      if (!exports.has(name)) throw new Error(`engine.wasm does not export ${name}`);
    }
    const [mathImport] = WebAssembly.Module.imports(mathlibModule);
    if (mathImport) {
      throw new Error(
        `mathlib.wasm imports ${mathImport.module}.${mathImport.name}, it imports nothing`
      );
    }
  }
  function validateModelBytes(source) {
    const bytes = source instanceof ArrayBuffer ? new Uint8Array(source) : ArrayBuffer.isView(source) ? new Uint8Array(source.buffer, source.byteOffset, source.byteLength) : null;
    if (!bytes || bytes.length < 8) throw new Error("model.wasm is not a WASM binary");
    let pos = 8;
    const leb = (end) => {
      let value = 0;
      for (let shift = 0; shift < 35 && pos < end; shift += 7) {
        const byte = bytes[pos++];
        value += (byte & 127) * 2 ** shift;
        if (!(byte & 128)) return value;
      }
      throw new Error("invalid model.wasm section encoding");
    };
    while (pos < bytes.length) {
      const id = bytes[pos++];
      const length = leb(bytes.length);
      const end = pos + length;
      if (end > bytes.length) throw new Error("invalid model.wasm section length");
      if (id === 1) {
        const count = leb(end);
        for (let i = 0; i < count; i++) {
          if (bytes[pos++] !== 96) throw new Error("model.wasm has an unsupported type");
          for (let vector = 0; vector < 2; vector++) {
            const n = leb(end);
            for (let j = 0; j < n; j++) {
              if (pos >= end || ![127, 126, 125, 124, 123].includes(bytes[pos++])) {
                throw new Error("model.wasm has an unsupported reference type");
              }
            }
          }
        }
      }
      if (id === 13) throw new Error("model.wasm must not define exception tags");
      if ([4, 5, 8, 9, 11, 12].includes(id)) {
        if (id === 8 || id === 12 || leb(end) !== 0) {
          throw new Error("model.wasm must not define tables, memories, a start function, elements or data");
        }
      }
      if (id === 6) {
        const count = leb(end);
        for (let i = 0; i < count; i++) {
          const type = bytes[pos++];
          if (![127, 126, 125, 124].includes(type)) throw new Error("model.wasm has an unsupported global type");
          if (bytes[pos++] !== 0) throw new Error("model.wasm must not define mutable globals");
          const op = bytes[pos++];
          if (op === 65 || op === 66) {
            let n = 0;
            while (pos < end && n++ < 10 && bytes[pos++] & 128) {
            }
          } else if (op === 67) pos += 4;
          else if (op === 68) pos += 8;
          else throw new Error("model.wasm has an unsupported global initializer");
          if (pos >= end || bytes[pos++] !== 11) throw new Error("model.wasm has an invalid global initializer");
        }
      }
      pos = end;
    }
  }
  function checkModel(modelModule, mathlibModule) {
    const math = new Set(
      WebAssembly.Module.exports(mathlibModule).filter((e) => e.kind === "function").map((e) => e.name)
    );
    let memory = false;
    for (const { module, name, kind } of WebAssembly.Module.imports(modelModule)) {
      if (module === "env" && name === "memory" && kind === "memory") memory = true;
      else if (module === "math" && kind === "function" && math.has(name)) continue;
      else {
        throw new Error(
          `model.wasm imports ${module}.${name}, which is neither env.memory nor a function of mathlib`
        );
      }
    }
    if (!memory) throw new Error("model.wasm does not import env.memory");
  }
  function jsonText(value) {
    return typeof value === "string" ? value : JSON.stringify(value);
  }
  async function createRuntime({
    engine,
    mathlib,
    console: out = globalThis.console,
    maxMemoryPages = MAX_MEMORY_PAGES
  }) {
    const engineModule = await toModule(engine, "engine.wasm");
    const mathlibModule = await toModule(mathlib, "mathlib.wasm");
    checkRuntime(engineModule, mathlibModule);
    const minimum = memoryMinimum(engineModule);
    if (!Number.isSafeInteger(maxMemoryPages) || maxMemoryPages < minimum || maxMemoryPages > MAX_MEMORY_PAGES) {
      throw new Error(`maxMemoryPages must be between ${minimum} and ${MAX_MEMORY_PAGES}`);
    }
    return {
      /** The compiled modules, which a worker gets to create a runtime of its own. */
      modules: { engine: engineModule, mathlib: mathlibModule },
      /**
       * Instantiates a model (a `WebAssembly.Module` or the bytes of `model.wasm`) with its
       * description (the text of `model.json` or its object) and loads it into a new engine.
       * Asynchronous (a Promise of the instance): bytes compile with `WebAssembly.compile`, which
       * browsers require for modules above a few KB on the main thread.
       */
      async instantiate(modelWasm, description) {
        const descText = jsonText(description);
        const desc = typeof description === "string" ? JSON.parse(description) : description;
        if (desc.contract !== MODEL_CONTRACT) {
          throw new Error(`the description is of the contract ${desc.contract}, not ${MODEL_CONTRACT}`);
        }
        if (!(modelWasm instanceof WebAssembly.Module)) validateModelBytes(modelWasm);
        const modelModule = await toModule(modelWasm, "model.wasm");
        const history = WebAssembly.Module.customSections(modelModule, "sbmlsim.history");
        const nDelays = desc.n_delays ?? 0;
        if (history.length > 1 || history.length === 0 && nDelays !== 0) {
          throw new Error("model.wasm history sites differ from its description");
        }
        if (history.length === 1) {
          const info = new DataView(history[0]);
          if (info.byteLength !== 8 || info.getUint32(0, true) !== 1 || info.getUint32(4, true) !== nDelays) {
            throw new Error("model.wasm history metadata differs from its description");
          }
        }
        checkModel(modelModule, mathlibModule);
        const memory = new WebAssembly.Memory({ initial: minimum, maximum: maxMemoryPages });
        const wasi = wasiImports(engineModule, () => memory, out);
        const math = await WebAssembly.instantiate(mathlibModule, {});
        const model = await WebAssembly.instantiate(modelModule, {
          env: { memory },
          math: math.exports
        });
        let engineInstance;
        try {
          engineInstance = await WebAssembly.instantiate(engineModule, {
            env: { memory },
            model: model.exports,
            [WASI]: wasi.imports
          });
        } catch (error) {
          if (error instanceof WebAssembly.LinkError) {
            throw new Error(`engine.wasm does not instantiate with the model: ${error.message}`);
          }
          throw error;
        }
        const instance = new Instance(engineInstance.exports, memory, wasi.calls);
        if (engineInstance.exports._initialize) {
          instance.guard(() => engineInstance.exports._initialize());
        }
        instance.initializeCalls = wasi.calls();
        instance.load(descText, desc.id);
        return instance;
      }
    };
  }
  var preparedStates = /* @__PURE__ */ new WeakMap();
  var PreparedSimulation = class {
    constructor(token, state) {
      if (token !== preparedStates) throw new TypeError("use instance.prepare()");
      preparedStates.set(this, state);
    }
    get nValues() {
      return preparedStates.get(this).nValues;
    }
    close() {
      const state = preparedStates.get(this);
      if (state.closed) return;
      state.closed = true;
      state.instance.plans.delete(this);
      if (!state.instance.trapped) state.instance.guard(() => state.instance.exports.engine_release_prepared(state.handle));
    }
  };
  var Instance = class {
    constructor(exports, memory, calls) {
      this.exports = exports;
      this.memory = memory;
      this.calls = calls;
      this.trapped = false;
      this.encoder = new TextEncoder();
      this.decoder = new TextDecoder();
      this.inputPtr = 0;
      this.inputCapacity = 0;
      this.plans = /* @__PURE__ */ new Set();
    }
    /** The number of calls of the engine into WASI since the instance was created. */
    wasiCalls() {
      return this.calls();
    }
    /** The size of the memory of the instance in bytes. */
    memoryBytes() {
      return this.memory.buffer.byteLength;
    }
    /**
     * The counters of the integrator of the last simulation: an object of the names of
     * `STATS_NAMES` to numbers (`engine_stats_ptr`).
     */
    stats() {
      if (this.trapped) {
        throw new EngineError("internal", "the instance trapped before and is not used again");
      }
      const ptr = this.guard(() => this.exports.engine_stats_ptr()) >>> 0;
      const view = new DataView(this.memory.buffer, ptr, 8 * STATS_NAMES.length);
      return Object.fromEntries(STATS_NAMES.map((name, i) => [name, Number(view.getBigUint64(8 * i, true))]));
    }
    /** The address of the output of the last call in the memory (aligned to 8). */
    outputAddress() {
      return this.exports.engine_output_ptr() >>> 0;
    }
    /**
     * Simulates a request (an object of the contract `simulation/1` or its JSON text); gives
     * `{ columns, nRows, values }`, the values a `Float64Array` of the rows, row major. The result
     * does not depend on the simulations before. An error of the engine is an `EngineError`.
     */
    simulate(request) {
      const code = this.call(this.exports.engine_simulate, jsonText(request));
      const output = this.output();
      if (code === 0 && output.table) return output.table;
      if (code === 1 && output.error) throw output.error;
      throw new EngineError("internal", `engine_simulate returned ${code} with the output ${JSON.stringify(output.header)}`);
    }
    /** Prepares fixed request structure once; varying IDs define the raw value order. */
    prepare(request, { varying = [] } = {}) {
      if (!["engine_prepare", "engine_simulate_prepared", "engine_release_prepared"].every((name) => typeof this.exports[name] === "function")) {
        throw new UnsupportedError("this engine does not support prepared requests");
      }
      if (this.trapped) throw new EngineError("internal", "the instance trapped before and is not used again");
      if (this.plans.size >= 32) throw new EngineError("request", "prepared request limit reached (32)");
      const code = this.call(this.exports.engine_prepare, jsonText({ request: typeof request === "string" ? JSON.parse(request) : request, varying }));
      const output = this.output();
      if (code !== 0) throw output.error ?? new EngineError("internal", "engine_prepare failed without an error");
      const { prepared: handle, n_values: nValues } = output.header;
      if (!Number.isInteger(handle) || !Number.isSafeInteger(nValues) || nValues < 0) throw new EngineError("internal", "invalid prepared request output");
      const plan = new PreparedSimulation(preparedStates, { instance: this, handle, nValues, closed: false });
      this.plans.add(plan);
      return plan;
    }
    simulatePrepared(plan, values = []) {
      const state = preparedStates.get(plan);
      if (!state || state.instance !== this || state.closed) throw new EngineError("request", "foreign or closed prepared request");
      if (this.trapped) throw new EngineError("internal", "the instance trapped before and is not used again");
      if (!Array.isArray(values) && !(values instanceof Float64Array) || values.length !== state.nValues || values.findIndex((value) => !Number.isFinite(value)) >= 0) {
        throw new EngineError("request", `prepared request requires ${state.nValues} finite values`);
      }
      const ptr = this.inputBuffer(values.length * 8);
      const view = new DataView(this.memory.buffer, ptr, values.length * 8);
      for (let i = 0; i < values.length; i++) view.setFloat64(i * 8, values[i], true);
      const code = this.guard(() => this.exports.engine_simulate_prepared(state.handle, ptr, values.length));
      const output = this.output();
      if (code === 0 && output.table) return output.table;
      if (output.error) throw output.error;
      throw new EngineError("internal", "engine_simulate_prepared returned invalid output");
    }
    inputBuffer(length) {
      if (length > this.inputCapacity) {
        const capacity = Math.max(length, this.inputCapacity * 2, 256);
        const ptr = this.guard(() => this.exports.engine_alloc(capacity)) >>> 0;
        if (!ptr) throw new EngineError("internal", "engine_alloc found no memory for the input");
        if (this.inputCapacity) this.guard(() => this.exports.engine_free(this.inputPtr, this.inputCapacity));
        this.inputPtr = ptr;
        this.inputCapacity = capacity;
      }
      return this.inputPtr;
    }
    /** Loads the description (`engine_load`). */
    load(descText, id) {
      const before = this.calls();
      const code = this.call(this.exports.engine_load, descText);
      this.loadCalls = this.calls() - before;
      if (code !== 0) {
        const { error } = this.output();
        const message = error?.kind === "load" ? error.detail : `engine_load failed: ${error?.message}`;
        throw new Error(`the engine refuses the description of model ${id}: ${message}`);
      }
    }
    /** Calls `engine_load` or `engine_simulate` with `text` in a buffer of the engine. */
    call(fn, text) {
      if (this.trapped) {
        throw new EngineError("internal", "the instance trapped before and is not used again");
      }
      const bytes = this.encoder.encode(text);
      const ptr = this.inputBuffer(bytes.length);
      new Uint8Array(this.memory.buffer, ptr, bytes.length).set(bytes);
      const code = this.guard(() => fn(ptr, bytes.length));
      return code;
    }
    /**
     * Runs a call of the engine; whatever it throws marks the instance and is an `EngineError`
     * `internal` ([`trapError`]), as a trap of the wasmtime host, but `WasiExit`.
     */
    guard(call) {
      try {
        return call();
      } catch (error) {
        this.trapped = true;
        if (error instanceof WasiExit) throw error;
        let panic = null;
        try {
          const { error: e } = this.output();
          if (e?.engineKind === "panic") panic = e.detail;
        } catch {
        }
        throw trapError(error, panic);
      }
    }
    /** The output of the last call: `{ header, table }` or `{ header, error }`. */
    output() {
      const ptr = this.exports.engine_output_ptr() >>> 0;
      const len = this.exports.engine_output_len() >>> 0;
      const buffer = this.memory.buffer;
      if (len < 4 || ptr + len > buffer.byteLength) {
        throw new EngineError("internal", "the output of the engine is outside of the memory");
      }
      const view = new DataView(buffer, ptr, len);
      const headerLength = view.getUint32(0, true);
      if (4 + headerLength > len) {
        throw new EngineError("internal", "the header of the output of the engine is cut off");
      }
      const header = JSON.parse(this.decoder.decode(new Uint8Array(buffer, ptr + 4, headerLength)));
      if (header.error) {
        const { kind, message, flag, time } = header.error;
        const known = ["load", "request", "unsupported", "integration"];
        const error = known.includes(kind) ? new EngineError(kind, message, { flag, time }) : new EngineError("internal", kind === "panic" ? `engine.wasm panicked: ${message}` : `an error of the unknown kind ${kind}: ${message}`);
        error.detail = message;
        error.engineKind = kind;
        return { header, error };
      }
      if (header.prepared !== void 0) return { header };
      const { columns, n_rows: nRows } = header;
      const start = 4 + headerLength + 7 & ~7;
      const n = nRows * columns.length;
      if (start + 8 * n !== len) {
        throw new EngineError(
          "internal",
          `the output of the engine does not hold ${nRows} rows of ${columns.length} columns`
        );
      }
      let values;
      if (LITTLE_ENDIAN) {
        values = new Float64Array(buffer, ptr + start, n).slice();
      } else {
        values = new Float64Array(n);
        for (let i = 0; i < n; i += 1) values[i] = view.getFloat64(start + 8 * i, true);
      }
      return { header, table: { columns, nRows, values } };
    }
  };
  function trapError(error, panic = null) {
    const what = error instanceof WebAssembly.RuntimeError ? "trapped" : `failed (${error?.name ?? "Error"})`;
    const message = `engine.wasm ${what}: ${error?.message ?? String(error)}${panic ? ` (${panic})` : ""}`;
    const mapped = new EngineError("internal", message);
    mapped.cause = error;
    return mapped;
  }

  // pool.mjs
  var poolPlanStates = /* @__PURE__ */ new WeakMap();
  var preparedJob = /* @__PURE__ */ Symbol("prepared job");
  var PreparedPoolSimulation = class {
    constructor(token, state) {
      if (token !== poolPlanStates) throw new TypeError("use pool.prepare()");
      poolPlanStates.set(this, state);
    }
    get nValues() {
      return poolPlanStates.get(this).nValues;
    }
    close() {
      poolPlanStates.get(this).close();
    }
  };
  var WORKER_URL = new URL("./worker.mjs", "https://sbmlsim-core.invalid/");
  function defaultWorkers() {
    const cores = globalThis.navigator?.hardwareConcurrency ?? 2;
    return Math.max(1, Math.min(8, cores - 1));
  }
  var trampolines = /* @__PURE__ */ new Map();
  function browserWorker(url) {
    if (typeof globalThis.Worker !== "function") {
      throw new Error("this host has no Web Workers: give the pool a createWorker");
    }
    const script = new URL(url, globalThis.location?.href);
    const origin = globalThis.location?.origin;
    if (origin === void 0 || origin === "null" || script.origin === origin) {
      return new globalThis.Worker(script, { type: "module" });
    }
    let trampoline = trampolines.get(script.href);
    if (trampoline === void 0) {
      const source = `import ${JSON.stringify(script.href)};
`;
      trampoline = URL.createObjectURL(new Blob([source], { type: "text/javascript" }));
      trampolines.set(script.href, trampoline);
    }
    return new globalThis.Worker(trampoline, { type: "module" });
  }
  function listen(worker, onMessage, onError) {
    if (typeof worker.on === "function") {
      worker.on("message", onMessage);
      worker.on("error", onError);
      worker.on("messageerror", onError);
      worker.on("exit", (code) => onError(new Error(`the worker exited (${code})`)));
    } else {
      worker.addEventListener("message", (event) => onMessage(event.data));
      worker.addEventListener("error", (event) => {
        event.preventDefault?.();
        onError(event.error ?? new Error(event.message ?? "the worker failed"));
      });
      worker.addEventListener("messageerror", () => onError(new Error("a message of the worker could not be read")));
    }
  }
  function remoteError(fields) {
    if (fields.kind !== void 0 && ["EngineError", "UnsupportedError"].includes(fields.name)) {
      const error2 = new EngineError(fields.kind, fields.detail ?? fields.message, { flag: fields.flag, time: fields.time });
      if (fields.message) error2.message = fields.message;
      error2.name = fields.name;
      return error2;
    }
    const error = new Error(fields.message);
    if (fields.name) error.name = fields.name;
    return error;
  }
  async function createPool({
    engine,
    mathlib,
    workers = defaultWorkers(),
    createWorker = browserWorker,
    workerUrl = WORKER_URL,
    maxWorkers = 32,
    maxQueuedJobs = 1024,
    maxQueuedBytes = 64 * 2 ** 20,
    maxModelBytes = 128 * 2 ** 20,
    maxModels = 128,
    maxWorkerModels = 8,
    maxMemoryPages,
    maxTotalMemoryBytes = 2 ** 31,
    maxModuleBytes = 32 * 2 ** 20,
    maxDescriptionBytes = 16 * 2 ** 20,
    maxRequestBytes = 16 * 2 ** 20,
    deadlineMs = 3e4,
    hardAbort = false
  }) {
    if (!Number.isInteger(workers) || workers < 1) {
      throw new Error(`a pool needs at least one worker, not ${workers}`);
    }
    for (const [name, value] of Object.entries({ maxWorkers, maxQueuedJobs, maxQueuedBytes, maxModelBytes, maxModels, maxWorkerModels, maxTotalMemoryBytes, maxModuleBytes, maxDescriptionBytes, maxRequestBytes })) {
      if (!Number.isSafeInteger(value) || value < 1) throw new Error(`${name} must be a positive safe integer`);
    }
    if (workers > maxWorkers) throw new Error(`workers exceeds maxWorkers (${maxWorkers})`);
    const checkDeadline = (value) => {
      if (!Number.isFinite(value) || value < 0 || value > 2147483647) throw new Error("deadlineMs must be between 0 and 2147483647");
    };
    checkDeadline(deadlineMs);
    const checkModuleSize = (module) => {
      if (!(module instanceof WebAssembly.Module) && module?.byteLength > maxModuleBytes) {
        throw new Error(`module exceeds maxModuleBytes (${maxModuleBytes})`);
      }
    };
    checkModuleSize(engine);
    checkModuleSize(mathlib);
    const runtime = await createRuntime({ engine, mathlib, maxMemoryPages });
    const retained = workers * Math.min(maxWorkerModels, maxModels);
    const aggregatePages = Math.floor(maxTotalMemoryBytes / 65536 / retained);
    const sections = WebAssembly.Module.customSections(runtime.modules.engine, ENGINE_SECTION);
    const minimum = JSON.parse(new TextDecoder().decode(sections[0])).memory_minimum_pages;
    if (aggregatePages < minimum) {
      throw new Error(`maxTotalMemoryBytes (${maxTotalMemoryBytes}) cannot fit ${retained} retained instances at the engine minimum (${minimum} pages)`);
    }
    if (maxMemoryPages !== void 0 && maxMemoryPages > aggregatePages) {
      throw new Error(`maxMemoryPages (${maxMemoryPages}) exceeds maxTotalMemoryBytes (${maxTotalMemoryBytes}) for ${retained} retained instances; lower workers/maxWorkerModels or raise the total budget`);
    }
    maxMemoryPages ??= Math.min(MAX_MEMORY_PAGES, aggregatePages);
    const checkTextSize = (text, maximum, name) => {
      if (typeof text !== "string") throw new Error(`${name} must serialize to JSON text`);
      let bytes = 0;
      for (let i = 0; i < text.length; i++) {
        const code = text.charCodeAt(i);
        if (code < 128) bytes += 1;
        else if (code < 2048) bytes += 2;
        else if (code >= 55296 && code <= 56319 && text.charCodeAt(i + 1) >= 56320 && text.charCodeAt(i + 1) <= 57343) {
          bytes += 4;
          i++;
        } else bytes += 3;
        if (bytes > maximum) throw new Error(`${name} exceeds its byte limit (${maximum})`);
      }
      return bytes;
    };
    const models = /* @__PURE__ */ new Map();
    const plans = /* @__PURE__ */ new Map();
    let nextPlanId = 1;
    const loadingModels = /* @__PURE__ */ new Set();
    const slots = [];
    const queue = /* @__PURE__ */ new Map();
    let queuedBytes = 0;
    let modelBytes = 0;
    const dequeue = (job) => {
      if (queue.delete(job.id)) queuedBytes -= job.bytes;
    };
    let pumping = false;
    let nextId = 1;
    let closed = false;
    const settle = (job, error, result, failed = true) => {
      if (job.settled) return;
      job.settled = true;
      job.signal?.removeEventListener("abort", job.onAbort);
      if (failed) job.reject(error);
      else job.resolve(result);
    };
    const stop = (worker) => {
      try {
        Promise.resolve(worker.terminate()).catch(() => {
        });
      } catch {
      }
    };
    const spawn = () => {
      const slot = { worker: createWorker(workerUrl), models: /* @__PURE__ */ new Set(), plans: /* @__PURE__ */ new Set(), job: null, dead: false };
      try {
        listen(
          slot.worker,
          (message) => {
            const job = slot.job;
            if (slot.dead || !job || message?.id !== job.id) return;
            clearTimeout(job.timer);
            slot.job = null;
            if (message.type === "prepared") {
              settle(job, null, message.nValues, false);
            } else if (message.type === "result") {
              settle(job, null, { columns: message.columns, nRows: message.nRows, values: message.values }, false);
            } else {
              settle(job, remoteError(message.error ?? {}));
            }
            pump();
          },
          (error) => retire(slot, error)
        );
        slot.worker.postMessage({ type: "init", ...runtime.modules, maxMemoryPages });
      } catch (error) {
        slot.dead = true;
        stop(slot.worker);
        throw error;
      }
      slots.push(slot);
      return slot;
    };
    const retire = (slot, error) => {
      if (slot.dead) return;
      slot.dead = true;
      const index = slots.indexOf(slot);
      if (index >= 0) slots.splice(index, 1);
      stop(slot.worker);
      if (slot.job) {
        clearTimeout(slot.job.timer);
        settle(slot.job, error);
      }
      slot.job = null;
      pump();
    };
    const pump = () => {
      if (pumping) return;
      pumping = true;
      try {
        while (!closed && queue.size > 0) {
          const job = queue.values().next().value;
          const key = job.key;
          const idle = slots.filter((slot2) => !slot2.job);
          let slot = idle.find((s) => s.models.has(key)) ?? idle[0];
          if (!slot && slots.length >= workers) return;
          dequeue(job);
          try {
            slot ??= spawn();
            slot.job = job;
            if (job.deadlineMs > 0) job.timer = setTimeout(() => {
              const error = new Error(`worker deadline exceeded (${job.deadlineMs} ms)`);
              error.name = "TimeoutError";
              retire(slot, error);
            }, job.deadlineMs);
            if (!slot.models.has(job.key)) {
              if (slot.models.size >= maxWorkerModels) {
                const oldest = slot.models.values().next().value;
                slot.worker.postMessage({ type: "removeModel", key: oldest });
                slot.models.delete(oldest);
                for (const id of slot.plans) if (plans.get(id)?.key === oldest) slot.plans.delete(id);
              }
              const { module, description } = models.get(job.key);
              slot.worker.postMessage({ type: "model", key: job.key, module, description });
            }
            slot.models.delete(job.key);
            slot.models.add(job.key);
            if (job.prepared) {
              const state = job.prepared;
              const definition = slot.plans.has(state.id) ? {} : { request: state.request, varying: state.varying };
              slot.worker.postMessage({ type: job.type, id: job.id, key: job.key, planId: state.id, values: job.values, ...definition }, job.values ? [job.values.buffer] : []);
              slot.plans.add(state.id);
            } else slot.worker.postMessage({ type: "simulate", id: job.id, key: job.key, request: job.request });
          } catch (error) {
            settle(job, error);
            if (slot) retire(slot, error);
          }
        }
      } finally {
        pumping = false;
      }
    };
    const api = {
      /** The most workers of the pool. */
      size: workers,
      /** Effective per-instance and aggregate linear-memory policies. */
      maxMemoryPages,
      maxTotalMemoryBytes,
      /** The compiled modules of the runtime. */
      modules: runtime.modules,
      /**
       * Adds a model under a key: `model` is a `WebAssembly.Module` or the bytes of `model.wasm`
       * (compiled here once), `description` the text of `model.json` or its object. A Promise;
       * a key which exists already is refused.
       */
      async addModel(key, model, description) {
        if (closed) throw new Error("the pool is closed");
        if (models.has(key) || loadingModels.has(key)) throw new Error(`the pool has a model ${key} already`);
        if (models.size + loadingModels.size >= maxModels) throw new Error(`model limit reached (${maxModels})`);
        checkModuleSize(model);
        if (!(model instanceof WebAssembly.Module)) validateModelBytes(model);
        const text = typeof description === "string" ? description : JSON.stringify(description);
        const bytes = checkTextSize(text, maxDescriptionBytes, "description") + (model instanceof WebAssembly.Module ? 0 : model.byteLength);
        if (modelBytes + bytes > maxModelBytes) throw new Error(`registered models exceed maxModelBytes (${maxModelBytes})`);
        loadingModels.add(key);
        modelBytes += bytes;
        let registered = false;
        try {
          const module = model instanceof WebAssembly.Module ? model : await WebAssembly.compile(model);
          if (closed) throw new Error("the pool is closed");
          models.set(key, { module, description: text, bytes });
          registered = true;
        } finally {
          if (loadingModels.delete(key) && !registered) modelBytes -= bytes;
        }
      },
      /** Removes an idle model; queued/running users must complete or abort first. */
      removeModel(key) {
        if (loadingModels.has(key) || [...queue.values()].some((job) => job.key === key) || slots.some((slot) => slot.job?.key === key)) {
          throw new Error(`model ${key} is in use`);
        }
        const model = models.get(key);
        if (!model) return false;
        for (const state of [...plans.values()]) if (state.key === key) state.close();
        models.delete(key);
        modelBytes -= model.bytes;
        for (const slot of [...slots]) {
          if (!slot.models.delete(key)) continue;
          try {
            slot.worker.postMessage({ type: "removeModel", key });
          } catch (error) {
            retire(slot, error);
          }
        }
        return true;
      },
      /** Whether the pool has a model under a key. */
      hasModel(key) {
        return models.has(key);
      },
      async prepare(key, request, { varying = [] } = {}) {
        if (closed) throw new Error("the pool is closed");
        if (!models.has(key)) throw new Error(`the pool has no model ${key}`);
        if (plans.size >= 32) throw new Error("prepared request limit reached (32)");
        const text = typeof request === "string" ? request : JSON.stringify(request);
        const bytes = checkTextSize(text, maxRequestBytes, "request");
        const varyingText = JSON.stringify(varying);
        const varyingBytes = checkTextSize(varyingText, maxRequestBytes, "varying");
        if (modelBytes + bytes + varyingBytes > maxModelBytes) throw new Error(`prepared requests exceed maxModelBytes (${maxModelBytes})`);
        const state = { owner: api, id: nextPlanId++, key, request: text, varying: JSON.parse(varyingText), bytes: bytes + varyingBytes, nValues: 0, closed: false };
        const plan = new PreparedPoolSimulation(poolPlanStates, state);
        state.close = () => {
          if (state.closed) return;
          state.closed = true;
          plans.delete(state.id);
          modelBytes -= state.bytes;
          const error = new Error("the prepared request is closed");
          for (const job of [...queue.values()]) if (job.prepared === state) {
            dequeue(job);
            settle(job, error);
          }
          for (const slot of [...slots]) if (slot.plans.delete(state.id)) {
            try {
              slot.worker.postMessage({ type: "releasePrepared", key, planId: state.id });
            } catch (error2) {
              retire(slot, error2);
            }
          }
        };
        plans.set(state.id, state);
        modelBytes += state.bytes;
        try {
          state.nValues = await api.simulate(key, { [preparedJob]: { type: "prepare", prepared: state } });
          if (state.closed || closed) throw new Error("the prepared request is closed");
          return plan;
        } catch (error) {
          state.close();
          throw error;
        }
      },
      simulatePrepared(plan, values = [], options = {}) {
        const state = poolPlanStates.get(plan);
        if (!state || state.owner !== api || state.closed) return Promise.reject(new Error("foreign or closed prepared request"));
        if (!Array.isArray(values) && !(values instanceof Float64Array) || values.length !== state.nValues || values.findIndex((value) => !Number.isFinite(value)) >= 0) {
          return Promise.reject(new Error(`prepared request requires ${state.nValues} finite values`));
        }
        return api.simulate(state.key, { [preparedJob]: { type: "simulatePrepared", prepared: state, values: Float64Array.from(values) } }, options);
      },
      /**
       * Simulates a request of `simulation/1` (an object or its JSON text) on the model `key` in
       * a worker; resolves to `{ columns, nRows, values }` as `instance.simulate` of the loader,
       * the values a `Float64Array`. `signal` (an `AbortSignal`) drops the request, see the module
       * comment.
       */
      simulate(key, request, { signal, deadlineMs: jobDeadline = deadlineMs, hardAbort: jobHardAbort = hardAbort } = {}) {
        return new Promise((resolve, reject) => {
          if (closed) throw new Error("the pool is closed");
          if (!models.has(key)) throw new Error(`the pool has no model ${key}`);
          checkDeadline(jobDeadline);
          if (queue.size >= maxQueuedJobs) throw new Error(`pool queue limit reached (${maxQueuedJobs})`);
          signal?.throwIfAborted();
          const prepared = request?.[preparedJob];
          const text = prepared ? "" : typeof request === "string" ? request : JSON.stringify(request, (_key, value) => {
            if (typeof value === "function" || typeof value === "symbol") {
              throw new DOMException("request values could not be cloned", "DataCloneError");
            }
            return value;
          });
          const bytes = prepared ? prepared.values?.byteLength ?? prepared.prepared.bytes : checkTextSize(text, maxRequestBytes, "request");
          if (bytes > maxRequestBytes) throw new Error(`request exceeds its byte limit (${maxRequestBytes})`);
          if (queuedBytes + bytes > maxQueuedBytes) throw new Error(`pool queue exceeds maxQueuedBytes (${maxQueuedBytes})`);
          const job = { ...prepared, bytes, id: nextId++, key, request: text, signal, resolve, reject, settled: false, deadlineMs: jobDeadline };
          job.onAbort = () => {
            dequeue(job);
            settle(job, signal.reason);
            const slot = slots.find((slot2) => slot2.job === job);
            if (slot && jobHardAbort) retire(slot, signal.reason);
          };
          signal?.addEventListener("abort", job.onAbort, { once: true });
          queue.set(job.id, job);
          queuedBytes += bytes;
          pump();
        });
      },
      /**
       * Terminates the workers; every queued and running request is rejected, and the pool takes
       * no more.
       */
      close() {
        if (closed) return;
        closed = true;
        for (const state of [...plans.values()]) state.close();
        const error = new Error("the pool is closed");
        for (const job of queue.values()) settle(job, error);
        queue.clear();
        queuedBytes = 0;
        models.clear();
        loadingModels.clear();
        modelBytes = 0;
        for (const slot of slots.splice(0)) {
          slot.dead = true;
          stop(slot.worker);
          if (slot.job) {
            clearTimeout(slot.job.timer);
            settle(slot.job, error);
          }
        }
      }
    };
    return api;
  }

  // panels.mjs
  var CORE_VERSION = "0.24.0";
  var ENGINE = `sbmlsim-core ${CORE_VERSION}`;
  var BUNDLE_CONTRACT = "panelbundle/1";
  var PANEL_CONTRACT = "panel/1";
  var MAX_RUN_BYTES = 64 * 1024 * 1024;
  var MAX_ID = 128;
  var DEFAULT_LIMITS = { knobs: 8, grid: 64 };
  var REQUEST_FIELDS = /* @__PURE__ */ new Set(["values", "grid", "variant", "draws"]);
  var HttpError = class extends Error {
    constructor(status, detail, headers = {}) {
      super(detail);
      this.status = status;
      this.headers = headers;
    }
  };
  function json(status, body, headers = {}) {
    return new Response(JSON.stringify(body), {
      status,
      headers: { "content-type": "application/json", ...headers }
    });
  }
  function pyList(ids) {
    return `[${ids.map((id) => `'${id}'`).join(", ")}]`;
  }
  function knobValues(knob, n = knob.steps) {
    if (knob.kind === "choice") return knob.choices.map((choice) => choice.value);
    if (!Number.isInteger(n) || n < 1) throw new Error(`a grid over the knob '${knob.id}' needs a point, not ${n}`);
    const low = knob.min;
    const high = knob.max;
    if (n === 1) return [low];
    const last = n - 1;
    const inner = [];
    for (let i = 1; i < last; i += 1) inner.push((low * (last - i) + high * i) / last);
    return [low, ...inner, high];
  }
  function canonicalJson(value) {
    if (value === null || typeof value !== "object") {
      return JSON.stringify(value ?? null).replace(
        /[\u0080-￿]/g,
        (c) => `\\u${c.charCodeAt(0).toString(16).padStart(4, "0")}`
      );
    }
    if (Array.isArray(value)) return `[${value.map((v) => canonicalJson(v)).join(",")}]`;
    const keys = Object.keys(value).filter((key) => value[key] !== void 0).sort();
    return `{${keys.map((key) => `${canonicalJson(key)}:${canonicalJson(value[key])}`).join(",")}}`;
  }
  var SHA256_K = new Uint32Array([
    1116352408,
    1899447441,
    3049323471,
    3921009573,
    961987163,
    1508970993,
    2453635748,
    2870763221,
    3624381080,
    310598401,
    607225278,
    1426881987,
    1925078388,
    2162078206,
    2614888103,
    3248222580,
    3835390401,
    4022224774,
    264347078,
    604807628,
    770255983,
    1249150122,
    1555081692,
    1996064986,
    2554220882,
    2821834349,
    2952996808,
    3210313671,
    3336571891,
    3584528711,
    113926993,
    338241895,
    666307205,
    773529912,
    1294757372,
    1396182291,
    1695183700,
    1986661051,
    2177026350,
    2456956037,
    2730485921,
    2820302411,
    3259730800,
    3345764771,
    3516065817,
    3600352804,
    4094571909,
    275423344,
    430227734,
    506948616,
    659060556,
    883997877,
    958139571,
    1322822218,
    1537002063,
    1747873779,
    1955562222,
    2024104815,
    2227730452,
    2361852424,
    2428436474,
    2756734187,
    3204031479,
    3329325298
  ]);
  function sha256Bytes(bytes) {
    const length = bytes.length;
    const padded = new Uint8Array(length + 9 + 63 >> 6 << 6);
    padded.set(bytes);
    padded[length] = 128;
    const view = new DataView(padded.buffer);
    view.setUint32(padded.length - 8, Math.floor(length / 536870912));
    view.setUint32(padded.length - 4, length << 3 >>> 0);
    const h = new Uint32Array([1779033703, 3144134277, 1013904242, 2773480762, 1359893119, 2600822924, 528734635, 1541459225]);
    const w = new Uint32Array(64);
    const rotr = (x, n) => x >>> n | x << 32 - n;
    for (let offset = 0; offset < padded.length; offset += 64) {
      for (let i = 0; i < 16; i += 1) w[i] = view.getUint32(offset + 4 * i);
      for (let i = 16; i < 64; i += 1) {
        const s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ w[i - 15] >>> 3;
        const s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ w[i - 2] >>> 10;
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
      }
      let [a, b, c, d, e, f, g, k] = h;
      for (let i = 0; i < 64; i += 1) {
        const t1 = k + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + (e & f ^ ~e & g) + SHA256_K[i] + w[i];
        const t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + (a & b ^ a & c ^ b & c);
        k = g;
        g = f;
        f = e;
        e = d + t1 >>> 0;
        d = c;
        c = b;
        b = a;
        a = t1 + t2 >>> 0;
      }
      h[0] += a;
      h[1] += b;
      h[2] += c;
      h[3] += d;
      h[4] += e;
      h[5] += f;
      h[6] += g;
      h[7] += k;
    }
    const digest = new Uint8Array(32);
    const out = new DataView(digest.buffer);
    h.forEach((word, i) => out.setUint32(4 * i, word));
    return digest;
  }
  async function sha256Hex(text) {
    const bytes = new TextEncoder().encode(text);
    const digest = globalThis.crypto?.subtle ? new Uint8Array(await globalThis.crypto.subtle.digest("SHA-256", bytes)) : sha256Bytes(bytes);
    return Array.from(digest, (byte) => byte.toString(16).padStart(2, "0")).join("");
  }
  function route(path) {
    const parts = path.split("/").slice(1).map((part) => decodeURIComponent(part));
    if (parts.length < 4 || parts[0] !== "api" || parts[1] !== "models") return null;
    const model = parts[2];
    if (parts.length === 4 && parts[3] === "app") return { name: "app", method: "GET", model };
    if (parts[3] !== "panels" || parts.length < 5) return null;
    const panel = parts[4];
    if (parts.length === 5) return { name: "panel", method: "GET", model, panel };
    if (parts.length === 6 && parts[5] === "data") return { name: "data", method: "GET", model, panel };
    if (parts.length === 6 && parts[5] === "simulate") return { name: "simulate", method: "POST", model, panel };
    return null;
  }
  async function readInput(input, init) {
    const request = typeof input === "object" && input !== null && typeof input.method === "string" && !(input instanceof URL) ? input : null;
    const url = String(request ? request.url : input);
    const method = String(init.method ?? request?.method ?? "GET").toUpperCase();
    const signal = init.signal ?? request?.signal ?? void 0;
    let body = init.body ?? null;
    if (body === null && request && method !== "GET" && method !== "HEAD") body = await request.text();
    else if (body !== null && typeof body !== "string") body = await new Response(body).text();
    return { url, method, signal, body };
  }
  function readBody(text) {
    let body;
    try {
      body = text === null || text === "" ? {} : JSON.parse(text);
    } catch {
      throw new HttpError(422, "The body of a simulation is no JSON.");
    }
    if (body === null || typeof body !== "object" || Array.isArray(body)) {
      throw new HttpError(422, "The body of a simulation is no JSON object.");
    }
    const unknown = Object.keys(body).filter((key) => !REQUEST_FIELDS.has(key));
    if (unknown.length > 0) {
      throw new HttpError(422, `The body of a simulation has the unknown fields ${pyList(unknown.sort())}; its fields are ['draws', 'grid', 'values', 'variant'].`);
    }
    const mapping = (name) => {
      const value = body[name] ?? {};
      if (typeof value !== "object" || Array.isArray(value)) throw new HttpError(422, `The field '${name}' is no mapping of knob ids.`);
      for (const id of Object.keys(value)) {
        if (id.length > MAX_ID) throw new HttpError(422, `A knob id of the field '${name}' is longer than ${MAX_ID} characters.`);
      }
      return value;
    };
    const values = mapping("values");
    for (const [id, value] of Object.entries(values)) {
      if (typeof value !== "number" || !Number.isFinite(value)) {
        throw new HttpError(422, `The value of the knob '${id}' is no finite number.`);
      }
    }
    const grid = mapping("grid");
    for (const [id, n] of Object.entries(grid)) {
      if (!Number.isInteger(n) || n < 1) throw new HttpError(422, `The points of the grid over the knob '${id}' are no whole number of at least 1.`);
    }
    const variant = body.variant ?? "";
    if (typeof variant !== "string" || variant.length > MAX_ID) {
      throw new HttpError(422, `The variant is no text of at most ${MAX_ID} characters.`);
    }
    const draws = body.draws ?? 0;
    if (!Number.isInteger(draws) || draws < 0) {
      throw new HttpError(422, `The draws of the uncertainty are no whole number of at least 0.`);
    }
    return { values, grid, variant, draws };
  }
  function checkLimits(values, grid, limits) {
    for (const [name, mapping] of [["values", values], ["grid", grid]]) {
      const n = Object.keys(mapping).length;
      if (n > limits.knobs) throw new HttpError(422, `The field '${name}' names ${n} knobs, over the limit of ${limits.knobs}.`);
    }
    let points = 1;
    for (const [id, n] of Object.entries(grid)) {
      if (n > limits.grid) throw new HttpError(422, `The grid of ${n} points over the knob '${id}' is over the limit of ${limits.grid} points of a request.`);
      points *= n;
    }
    if (points > limits.grid) {
      throw new HttpError(422, `The grid of ${points} points is over the limit of ${limits.grid} points of a request.`);
    }
  }
  function checkKnobs(panelId, knobs, values, grid) {
    const ids = knobs.map((knob) => knob.id);
    const unknown = [.../* @__PURE__ */ new Set([...Object.keys(values), ...Object.keys(grid)])].filter((id) => !ids.includes(id)).sort();
    if (unknown.length > 0) {
      throw new HttpError(422, `${pyList(unknown)} are no knobs of the panel '${panelId}', which has ${pyList([...ids].sort())}.`);
    }
    const byId = new Map(knobs.map((knob) => [knob.id, knob]));
    for (const [id, value] of Object.entries(values)) {
      const knob = byId.get(id);
      let taken = "";
      if (knob.kind === "choice") {
        const choices = knob.choices.map((choice) => choice.value);
        if (!choices.includes(value)) taken = `none of the choices [${choices.join(", ")}]`;
      } else if (knob.min != null && knob.max != null && !(knob.min <= value && value <= knob.max)) {
        taken = `outside the bounds [${knob.min}, ${knob.max}]`;
      }
      if (taken) throw new HttpError(422, `The value ${value} of the knob '${id}' of the panel '${panelId}' is ${taken}.`);
    }
    for (const [id, n] of Object.entries(grid)) {
      const knob = byId.get(id);
      if (knob.kind === "choice" && n !== knob.choices.length) {
        throw new HttpError(
          422,
          `A grid of ${n} points over the knob '${id}' of the panel '${panelId}' is none it runs: a grid over it runs its ${knob.choices.length} choices.`
        );
      }
    }
  }
  function abortable(promise, signal) {
    if (!signal) return promise;
    signal.throwIfAborted();
    return new Promise((resolve, reject) => {
      const onAbort = () => reject(signal.reason);
      signal.addEventListener("abort", onAbort, { once: true });
      promise.then(
        (value) => {
          signal.removeEventListener("abort", onAbort);
          resolve(value);
        },
        (error) => {
          signal.removeEventListener("abort", onAbort);
          reject(error);
        }
      );
    });
  }
  function product(axes) {
    let points = [[]];
    for (const axis of axes) {
      const next = [];
      for (const point of points) for (const value of axis.values) next.push([...point, value]);
      points = next;
    }
    return points;
  }
  var Runs = class {
    constructor(maxBytes, maxEntries, maxInflight) {
      this.maxEntries = maxEntries;
      this.maxInflight = maxInflight;
      this.maxBytes = maxBytes;
      this.bytes = 0;
      this.inflight = /* @__PURE__ */ new Map();
      this.done = /* @__PURE__ */ new Map();
    }
    /**
     * The results of the run of `key` (one per point, `null` for a point which failed): a stored
     * one, a run in flight, or a new one of `start(signal)`, whose signal aborts when every
     * request which waits for it aborted or when the run failed.
     */
    get(key, start, signal) {
      signal?.throwIfAborted();
      const stored = this.done.get(key);
      if (stored) {
        this.done.delete(key);
        this.done.set(key, stored);
        return Promise.resolve(stored.results);
      }
      let entry = this.inflight.get(key);
      if (!entry) {
        if (this.inflight.size >= this.maxInflight) return Promise.reject(new Error("panel run limit reached"));
        const controller = new AbortController();
        entry = { controller, waiters: 0 };
        entry.promise = start(controller.signal).then(
          (results) => {
            if (this.inflight.get(key) === entry) this.inflight.delete(key);
            this.store(key, results);
            return results;
          },
          (error) => {
            if (this.inflight.get(key) === entry) this.inflight.delete(key);
            controller.abort(error);
            throw error;
          }
        );
        entry.promise.catch(() => {
        });
        this.inflight.set(key, entry);
      }
      entry.waiters += 1;
      if (!signal) return entry.promise;
      return new Promise((resolve, reject) => {
        const onAbort = () => {
          entry.waiters -= 1;
          if (entry.waiters === 0) {
            if (this.inflight.get(key) === entry) this.inflight.delete(key);
            entry.controller.abort(signal.reason);
          }
          reject(signal.reason);
        };
        if (signal.aborted) {
          onAbort();
          return;
        }
        signal.addEventListener("abort", onAbort, { once: true });
        entry.promise.then(
          (results) => {
            signal.removeEventListener("abort", onAbort);
            resolve(results);
          },
          (error) => {
            signal.removeEventListener("abort", onAbort);
            reject(error);
          }
        );
      });
    }
    /** Keeps the results of a run, evicting the least recently used ones beyond the bound. */
    store(key, results) {
      if (results.some((result) => result === null)) return;
      const bytes = 256 + key.length * 2 + results.reduce((sum, result) => sum + 128 + result.values.byteLength + JSON.stringify(result.columns).length * 2, 0);
      if (bytes > this.maxBytes) return;
      const previous = this.done.get(key);
      if (previous) {
        this.bytes -= previous.bytes;
        this.done.delete(key);
      }
      this.done.set(key, { results, bytes });
      this.bytes += bytes;
      for (const [old, entry] of this.done) {
        if (this.bytes <= this.maxBytes && this.done.size <= this.maxEntries) break;
        this.done.delete(old);
        this.bytes -= entry.bytes;
      }
    }
  };
  var MAX_DOSES = 1e3;
  function doseTimes(interval, duration, end) {
    if (!(interval > 0)) return [];
    const until = duration === null || duration === void 0 ? end : Math.min(duration, end);
    const times = [];
    for (let k = 1; k <= MAX_DOSES; k += 1) {
      const time = k * interval;
      if (!(time < until)) break;
      times.push(time);
    }
    return times;
  }
  function doseChanges(run, state) {
    const changes = [];
    for (const dosing of run.dosing ?? []) {
      const dose = state[dosing.dose] * run.knobs[dosing.dose].factor;
      const interval = state[dosing.interval] * dosing.interval_factor;
      const duration = dosing.duration === null || dosing.duration === void 0 ? null : state[dosing.duration] * dosing.duration_factor;
      for (const time of doseTimes(interval, duration, run.end)) {
        changes.push({ time, values: [], increments: [[dosing.target, dose]] });
      }
    }
    return changes.sort((a, b) => a.time - b.time);
  }
  function outputRows(result) {
    if (result.columns[0] !== "time") throw new Error(`a result has the first column ${result.columns[0]}, not time`);
    const width = result.columns.length;
    const time = (r) => result.values[r * width];
    const rows = [];
    for (let r = 0; r < result.nRows; r += 1) {
      if (r + 1 < result.nRows && time(r + 1) === time(r)) continue;
      rows.push(r);
    }
    return rows;
  }
  function resolveDraws(panel, uncertainty, asked) {
    if (!uncertainty || asked === 0) return 0;
    const max = uncertainty.max_draws ?? (uncertainty.draws ?? []).length;
    if (asked > max) {
      throw new HttpError(422, `The panel '${panel}' carries ${max} draws of its uncertainty, ${asked} were asked for.`);
    }
    return asked;
  }
  function percentile(sample, share) {
    const sorted = sample.filter((one) => one !== null && one !== void 0 && Number.isFinite(one)).sort((a, b) => a - b);
    if (sorted.length === 0) return null;
    const at = share * (sorted.length - 1);
    const low = Math.floor(at);
    const high = Math.ceil(at);
    return sorted[low] + (at - low) * (sorted[high] - sorted[low]);
  }
  function createPanelEndpoint(options) {
    const {
      bundle,
      runtime = {
        engine: new URL("./runtime/engine.wasm", "https://sbmlsim-core.invalid/"),
        mathlib: new URL("./runtime/mathlib.wasm", "https://sbmlsim-core.invalid/")
      },
      workers = defaultWorkers(),
      createWorker,
      pool: givenPool = null,
      fetch: load = (url) => globalThis.fetch(url),
      maxRunBytes = MAX_RUN_BYTES,
      maxRunEntries = 128,
      maxInflightRuns = 32
    } = options;
    if (bundle === void 0 || bundle === null) throw new Error("createPanelEndpoint needs the URL of a bundle");
    const bundleUrl = new URL(String(bundle), globalThis.location?.href);
    for (const [name, value] of Object.entries({ maxRunBytes, maxRunEntries, maxInflightRuns })) {
      if (!Number.isSafeInteger(value) || value < 1) throw new Error(`${name} must be a positive safe integer`);
    }
    const runs = new Runs(maxRunBytes, maxRunEntries, maxInflightRuns);
    const models = /* @__PURE__ */ new Map();
    let loading = null;
    let ownPool = null;
    let closed = false;
    const bytesOf = async (url) => {
      const response = await load(url);
      if (!response.ok) throw new Error(`${url} answers ${response.status}`);
      return new Uint8Array(await response.arrayBuffer());
    };
    const runtimeModule = async (source) => source instanceof WebAssembly.Module || ArrayBuffer.isView(source) || source instanceof ArrayBuffer ? source : bytesOf(new URL(String(source), globalThis.location?.href));
    const loaded = () => {
      if (!loading) {
        loading = (async () => {
          const response = await load(bundleUrl);
          if (!response.ok) throw new Error(`${bundleUrl} answers ${response.status}`);
          const data = await response.json();
          if (data?.contract !== BUNDLE_CONTRACT) {
            throw new Error(`the bundle is of the contract ${data?.contract}, not ${BUNDLE_CONTRACT}`);
          }
          let pool = givenPool;
          if (!pool) {
            const [engine, mathlib] = await Promise.all([runtimeModule(runtime.engine), runtimeModule(runtime.mathlib)]);
            if (closed) throw new Error("the endpoint is closed");
            ownPool = await createPool({ engine, mathlib, workers, ...createWorker ? { createWorker } : {} });
            if (closed) ownPool.close();
            pool = ownPool;
          }
          return { bundle: data, pool };
        })();
        loading.catch(() => {
          loading = null;
        });
      }
      return loading;
    };
    const modelOf = (bundleData, pool, key) => {
      let model = models.get(key);
      if (!model) {
        model = (async () => {
          const entry = bundleData.models?.[key];
          if (!entry) throw new Error(`the bundle has no model ${key}`);
          const dir = String(entry.compiled).replace(/\/+$/, "");
          const [wasm, text] = await Promise.all([
            bytesOf(new URL(`${dir}/model.wasm`, bundleUrl)),
            bytesOf(new URL(`${dir}/model.json`, bundleUrl)).then((bytes) => new TextDecoder().decode(bytes))
          ]);
          const description = JSON.parse(text);
          const settings = { ...entry.settings };
          const tolerances = entry.absolute_tolerances ?? {};
          if (Object.keys(tolerances).length > 0) {
            const scalar = settings.absolute_tolerance?.scalar;
            settings.absolute_tolerance = {
              per_state: description.states.map((state) => {
                const value = tolerances[state.id] ?? scalar;
                if (typeof value !== "number") {
                  throw new Error(`the model ${key} has no absolute tolerance of its state ${state.id}`);
                }
                return value;
              })
            };
          }
          const poolKey = `${bundleUrl.href}#${key}`;
          if (!pool.hasModel(poolKey)) {
            await pool.addModel(poolKey, wasm, text).catch((error) => {
              if (!pool.hasModel(poolKey)) throw error;
            });
          }
          return { poolKey, selections: entry.selections, settings };
        })();
        model.catch(() => models.delete(key));
        models.set(key, model);
      }
      return model;
    };
    const answersOf = (bundleData, model, variant) => {
      if (variant.length > MAX_ID) throw new HttpError(422, `The variant is no text of at most ${MAX_ID} characters.`);
      const answers = bundleData.answers?.[variant];
      if (!answers) {
        const variants = Object.keys(bundleData.answers ?? {}).filter((name) => name !== "").sort();
        throw new HttpError(422, `The model '${model}' has no variant '${variant}'; its variants are ${pyList(variants)}.`);
      }
      return answers;
    };
    const panelOf = (bundleData, model, panel, variant) => {
      const answer = answersOf(bundleData, model, variant).panels?.[panel];
      if (!answer) {
        const within = variant === "" ? "" : ` in its variant '${variant}'`;
        throw new HttpError(404, `'${panel}' is no panel of the model '${model}'${within}.`);
      }
      return answer;
    };
    const simulate = async (bundleData, pool, model, panel, text, signal) => {
      const started = globalThis.performance?.now() ?? Date.now();
      const { values, grid, variant, draws } = readBody(text);
      const answer = panelOf(bundleData, model, panel, variant);
      const declaration = answer.declaration;
      checkLimits(values, grid, { ...DEFAULT_LIMITS, ...declaration.limits });
      const knobs = declaration.axes ?? [];
      checkKnobs(panel, knobs, values, grid);
      const run = answer.run;
      const asked = resolveDraws(panel, run.uncertainty, draws);
      const info = await abortable(modelOf(bundleData, pool, run.model), signal);
      const axes = knobs.filter((knob) => knob.id in grid).map((knob) => ({
        id: knob.id,
        label: knob.label,
        unit: run.knobs?.[knob.id]?.unit ?? null,
        values: knobValues(knob, grid[knob.id])
      }));
      const fixed = {};
      for (const knob of knobs) if (!(knob.id in grid)) fixed[knob.id] = values[knob.id] ?? knob.default;
      const points = product(axes).map((point) => {
        const state = { ...fixed };
        axes.forEach((axis, i) => {
          state[axis.id] = point[i];
        });
        const preinit = new Map(run.preinit ?? []);
        for (const knob of knobs) {
          const spec = run.knobs?.[knob.id];
          if (!spec) throw new Error(`the run of the panel ${panel} has no knob ${knob.id}`);
          if (spec.target !== null && spec.target !== void 0) preinit.set(spec.target, state[knob.id] * spec.factor);
          if (spec.parameter !== null && spec.parameter !== void 0) preinit.set(spec.parameter, state[knob.id]);
        }
        return {
          preinit: [...preinit].sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0),
          changes: doseChanges(run, state)
        };
      });
      const runner = (list) => async (runSignal) => {
        const requests = list.map(({ preinit, changes }) => ({
          contract: "simulation/1",
          start: run.start,
          end: run.end,
          output: { uniform: { steps: run.points } },
          preinit,
          changes,
          steady_state: null,
          selections: info.selections,
          settings: info.settings
        }));
        return Promise.all(
          requests.map(
            (request) => pool.simulate(info.poolKey, request, { signal: runSignal }).catch((error) => {
              if (error instanceof EngineError) return null;
              throw error;
            })
          )
        );
      };
      const rowsOf = (list) => {
        const out = {};
        const failedAt = [];
        for (const curve of declaration.curves ?? []) out[curve.id] = [];
        list.forEach((result, index) => {
          if (result === null) failedAt.push(index);
          const rows = result === null ? [] : outputRows(result);
          if (result !== null && rows.length !== run.x.length) {
            throw new Error(`the run of the panel ${panel} has ${rows.length} rows, its x ${run.x.length}`);
          }
          for (const curve of declaration.curves ?? []) {
            const spec = run.curves?.[curve.id];
            if (!spec) throw new Error(`the run of the panel ${panel} has no curve ${curve.id}`);
            if (result === null) {
              out[curve.id].push(run.x.map(() => null));
              continue;
            }
            const width = result.columns.length;
            const column = result.columns.indexOf(spec.selection);
            if (column < 0) throw new Error(`the model ${run.model} does not select ${spec.selection}`);
            const row = new Array(rows.length);
            rows.forEach((r, i) => {
              const value = result.values[r * width + column] * spec.factor;
              row[i] = Number.isFinite(value) ? value : null;
            });
            out[curve.id].push(row);
          }
        });
        return { curves: out, failed: failedAt };
      };
      const runKey = JSON.stringify([run.model, run.start, run.end, run.points, points]);
      const results = await runs.get(runKey, runner(points), signal);
      const { curves, failed } = rowsOf(results);
      const uncertaintyOf = async () => {
        if (asked === 0) return null;
        const vectors = (run.uncertainty.draws ?? []).slice(0, asked);
        const parameters = run.uncertainty.parameters ?? [];
        const drawPoints = [];
        for (const point of points) {
          for (const vector of vectors) {
            const preinit = new Map(point.preinit);
            parameters.forEach((id, k) => preinit.set(id, vector[k]));
            drawPoints.push({
              preinit: [...preinit].sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0),
              changes: point.changes
            });
          }
        }
        const drawKey = JSON.stringify([run.model, run.start, run.end, run.points, drawPoints]);
        const drawResults = await runs.get(drawKey, runner(drawPoints), signal);
        const drawn = rowsOf(drawResults);
        const level = run.uncertainty.level ?? 0.9;
        const low = (1 - level) / 2;
        const high = (1 + level) / 2;
        const band = {};
        for (const curve of declaration.curves ?? []) {
          band[curve.id] = points.map((_, p) => {
            const sample = vectors.map((__, d) => drawn.curves[curve.id][p * asked + d]);
            const ends = { low: new Array(run.x.length), high: new Array(run.x.length) };
            for (let t = 0; t < run.x.length; t += 1) {
              const column = sample.map((row) => row[t]);
              ends.low[t] = percentile(column, low);
              ends.high[t] = percentile(column, high);
            }
            return ends;
          });
        }
        return {
          draws: asked,
          fit: run.uncertainty.fit ?? "",
          level,
          band,
          ...points.length === 1 ? { curves: drawn.curves } : {},
          failed: drawn.failed
        };
      };
      const uncertainty = await uncertaintyOf();
      const resolvedValues = {};
      for (const knob of knobs) if (!(knob.id in grid)) resolvedValues[knob.id] = fixed[knob.id];
      const resolvedGrid = {};
      for (const knob of knobs) if (knob.id in grid) resolvedGrid[knob.id] = grid[knob.id];
      const modelEntry = bundleData.models[run.model];
      const specHash = await sha256Hex(
        canonicalJson({
          contract: PANEL_CONTRACT,
          model,
          panel,
          variant,
          declaration,
          model_sha256: run.model_sha256,
          values: resolvedValues,
          grid: resolvedGrid,
          draws: asked,
          engine: ENGINE,
          settings: { settings: modelEntry.settings, absolute_tolerances: modelEntry.absolute_tolerances ?? {} }
        })
      );
      const finished = globalThis.performance?.now() ?? Date.now();
      return {
        contract: PANEL_CONTRACT,
        spec_hash: specHash,
        x: run.x,
        dims: axes,
        curves,
        units: run.units,
        provenance: {
          model,
          panel,
          model_sha256: run.model_sha256,
          variant: run.variant,
          engine: ENGINE,
          sbmlsim: bundleData.sbmlsim,
          seconds: (finished - started) / 1e3,
          created: (/* @__PURE__ */ new Date()).toISOString()
        },
        status: { ok: failed.length === 0, points: results.length, failed },
        ...uncertainty === null ? {} : { uncertainty }
      };
    };
    async function endpoint(input, init = {}) {
      const { url, method, signal, body } = await readInput(input, init ?? {});
      signal?.throwIfAborted();
      try {
        const at = url.indexOf("/api/");
        const [path, query = ""] = (at < 0 ? "" : url.slice(at)).split("#")[0].split("?");
        let target;
        try {
          target = route(path);
        } catch {
          target = null;
        }
        if (!target) throw new HttpError(404, "Not Found");
        if (method !== target.method) throw new HttpError(405, "Method Not Allowed", { allow: target.method });
        const { bundle: data, pool } = await abortable(
          loaded().catch((error) => {
            throw new HttpError(500, `The bundle ${bundleUrl} could not be loaded: ${error.message}`);
          }),
          signal
        );
        if (target.model !== data.model) throw new HttpError(404, `'${target.model}' is no model of this service.`);
        if (target.name === "simulate") {
          let answer2;
          try {
            answer2 = await simulate(data, pool, target.model, target.panel, body, signal);
          } catch (error) {
            if (error instanceof HttpError || error === signal?.reason) throw error;
            throw new HttpError(500, `The panel '${target.panel}' of the model '${target.model}' did not run: ${error.message}`);
          }
          return json(200, answer2);
        }
        const variant = new URLSearchParams(query).get("variant") ?? "";
        if (target.name === "app") {
          const app = answersOf(data, target.model, variant).app;
          if (app === null || app === void 0) throw new HttpError(404, `The model '${target.model}' declares no app.`);
          return json(200, app);
        }
        const answer = panelOf(data, target.model, target.panel, variant);
        return json(200, target.name === "panel" ? answer.declaration : answer.data);
      } catch (error) {
        if (signal?.aborted && error === signal.reason) throw error;
        if (error instanceof HttpError) return json(error.status, { detail: error.message }, error.headers);
        return json(500, { detail: `The endpoint failed: ${error?.message ?? error}` });
      }
    }
    endpoint.close = () => {
      closed = true;
      ownPool?.close();
    };
    return endpoint;
  }

  // inline.mjs
  var WORKER_SOURCE = '(() => {\n  // loader.mjs\n  var MODEL_CONTRACT = "model/2";\n  var MAX_MEMORY_PAGES = 32768;\n  var ENGINE_SECTION = "sbmlsim.engine";\n  var ENGINE_OUTPUT_FORMAT = 2;\n  var WASI = "wasi_snapshot_preview1";\n  var ENGINE_EXPORTS = [\n    "engine_alloc",\n    "engine_free",\n    "engine_load",\n    "engine_simulate",\n    "engine_output_ptr",\n    "engine_output_len",\n    "engine_stats_ptr"\n  ];\n  var STATS_NAMES = [\n    "steps",\n    "rhs_evals",\n    "jacobian_evals",\n    "jacobian_rhs_evals",\n    "jacobian_fallbacks",\n    "linear_setups",\n    "error_test_failures",\n    "nonlinear_iterations",\n    "nonlinear_convergence_failures",\n    "root_evals",\n    "integrations",\n    "probe_rhs_evals",\n    "steady_state_rhs_evals",\n    "linear_iterations",\n    "linear_convergence_failures",\n    "jacobian_vector_evals",\n    "preconditioner_evals",\n    "preconditioner_solves",\n    "sensitivity_rhs_evals",\n    "sensitivity_model_rhs_evals",\n    "sensitivity_jac_p_evals"\n  ];\n  var LITTLE_ENDIAN = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;\n  var EngineError = class extends Error {\n    constructor(kind, message, { flag, time } = {}) {\n      super(`${kind}: ${message}`);\n      this.name = "EngineError";\n      this.kind = kind;\n      this.detail = message;\n      if (flag !== void 0) this.flag = flag;\n      if (time !== void 0) this.time = time === null ? NaN : time;\n    }\n  };\n  var UnsupportedError = class extends EngineError {\n    constructor(message) {\n      super("unsupported", message);\n      this.name = "UnsupportedError";\n    }\n  };\n  var WasiExit = class extends Error {\n    constructor(code) {\n      super(`engine.wasm called proc_exit(${code})`);\n      this.name = "WasiExit";\n      this.code = code;\n    }\n  };\n  async function toModule(source, what) {\n    if (source instanceof WebAssembly.Module) return source;\n    try {\n      return await WebAssembly.compile(source);\n    } catch (error) {\n      throw new Error(`${what} does not compile: ${error.message}`);\n    }\n  }\n  function memoryMinimum(engineModule) {\n    const sections = WebAssembly.Module.customSections(engineModule, ENGINE_SECTION);\n    if (sections.length !== 1) {\n      throw new Error(\n        `engine.wasm has no custom section ${ENGINE_SECTION} with the minimum of its memory and the version of its output (the build script of sbmlsim-host adds it)`\n      );\n    }\n    const { memory_minimum_pages: pages, output_format: format } = JSON.parse(\n      new TextDecoder().decode(sections[0])\n    );\n    if (format === void 0) {\n      throw new Error(\n        `engine.wasm declares no output format in its section ${ENGINE_SECTION}, this loader reads ${ENGINE_OUTPUT_FORMAT}`\n      );\n    }\n    if (format !== ENGINE_OUTPUT_FORMAT) {\n      throw new Error(\n        `engine.wasm writes the output format ${format}, this loader reads ${ENGINE_OUTPUT_FORMAT}`\n      );\n    }\n    if (!Number.isInteger(pages) || pages < 0 || pages > MAX_MEMORY_PAGES) {\n      throw new Error(`engine.wasm declares a memory of ${pages} pages`);\n    }\n    return pages;\n  }\n  function wasiImports(engineModule, getMemory, out) {\n    let calls = 0;\n    const decoders = /* @__PURE__ */ new Map();\n    const lines = /* @__PURE__ */ new Map();\n    const view = () => new DataView(getMemory().buffer);\n    const write = (fd, bytes) => {\n      if (!decoders.has(fd)) decoders.set(fd, new TextDecoder());\n      let text = (lines.get(fd) ?? "") + decoders.get(fd).decode(bytes, { stream: true });\n      let newline;\n      while ((newline = text.indexOf("\\n")) >= 0) {\n        (fd === 2 ? out.error : out.log)(text.slice(0, newline));\n        text = text.slice(newline + 1);\n      }\n      lines.set(fd, text);\n    };\n    const provided = {\n      fd_write(fd, iovs, iovsLen, nwritten) {\n        calls += 1;\n        if (fd !== 1 && fd !== 2) return 8;\n        const memory = view();\n        let total = 0;\n        for (let i = 0; i < iovsLen; i += 1) {\n          const ptr = memory.getUint32(iovs + 8 * i, true);\n          const len = memory.getUint32(iovs + 8 * i + 4, true);\n          write(fd, new Uint8Array(getMemory().buffer, ptr, len));\n          total += len;\n        }\n        memory.setUint32(nwritten, total, true);\n        return 0;\n      },\n      proc_exit(code) {\n        calls += 1;\n        throw new WasiExit(code);\n      },\n      environ_get() {\n        calls += 1;\n        return 0;\n      },\n      environ_sizes_get(count, size) {\n        calls += 1;\n        const memory = view();\n        memory.setUint32(count, 0, true);\n        memory.setUint32(size, 0, true);\n        return 0;\n      },\n      clock_time_get() {\n        calls += 1;\n        throw new Error("engine.wasm called WASI clock_time_get: the engine reads no clock");\n      },\n      random_get() {\n        calls += 1;\n        throw new Error("engine.wasm called WASI random_get: the engine draws no random bytes");\n      }\n    };\n    const imports = {};\n    for (const { module, name, kind } of WebAssembly.Module.imports(engineModule)) {\n      if (module !== WASI || kind !== "function") continue;\n      imports[name] = provided[name] ?? (() => {\n        calls += 1;\n        throw new Error(`engine.wasm called WASI ${name}, which the loader does not provide`);\n      });\n    }\n    return { imports, calls: () => calls };\n  }\n  function checkRuntime(engineModule, mathlibModule) {\n    let memory = false;\n    for (const { module, name, kind } of WebAssembly.Module.imports(engineModule)) {\n      if (module === "env" && name === "memory" && kind === "memory") memory = true;\n      else if ((module === "model" || module === WASI) && kind === "function") continue;\n      else throw new Error(`engine.wasm imports ${module}.${name}`);\n    }\n    if (!memory) throw new Error("engine.wasm imports no memory");\n    const exports = new Set(WebAssembly.Module.exports(engineModule).map((e) => e.name));\n    for (const name of ENGINE_EXPORTS) {\n      if (!exports.has(name)) throw new Error(`engine.wasm does not export ${name}`);\n    }\n    const [mathImport] = WebAssembly.Module.imports(mathlibModule);\n    if (mathImport) {\n      throw new Error(\n        `mathlib.wasm imports ${mathImport.module}.${mathImport.name}, it imports nothing`\n      );\n    }\n  }\n  function validateModelBytes(source) {\n    const bytes = source instanceof ArrayBuffer ? new Uint8Array(source) : ArrayBuffer.isView(source) ? new Uint8Array(source.buffer, source.byteOffset, source.byteLength) : null;\n    if (!bytes || bytes.length < 8) throw new Error("model.wasm is not a WASM binary");\n    let pos = 8;\n    const leb = (end) => {\n      let value = 0;\n      for (let shift = 0; shift < 35 && pos < end; shift += 7) {\n        const byte = bytes[pos++];\n        value += (byte & 127) * 2 ** shift;\n        if (!(byte & 128)) return value;\n      }\n      throw new Error("invalid model.wasm section encoding");\n    };\n    while (pos < bytes.length) {\n      const id = bytes[pos++];\n      const length = leb(bytes.length);\n      const end = pos + length;\n      if (end > bytes.length) throw new Error("invalid model.wasm section length");\n      if (id === 1) {\n        const count = leb(end);\n        for (let i = 0; i < count; i++) {\n          if (bytes[pos++] !== 96) throw new Error("model.wasm has an unsupported type");\n          for (let vector = 0; vector < 2; vector++) {\n            const n = leb(end);\n            for (let j = 0; j < n; j++) {\n              if (pos >= end || ![127, 126, 125, 124, 123].includes(bytes[pos++])) {\n                throw new Error("model.wasm has an unsupported reference type");\n              }\n            }\n          }\n        }\n      }\n      if (id === 13) throw new Error("model.wasm must not define exception tags");\n      if ([4, 5, 8, 9, 11, 12].includes(id)) {\n        if (id === 8 || id === 12 || leb(end) !== 0) {\n          throw new Error("model.wasm must not define tables, memories, a start function, elements or data");\n        }\n      }\n      if (id === 6) {\n        const count = leb(end);\n        for (let i = 0; i < count; i++) {\n          const type = bytes[pos++];\n          if (![127, 126, 125, 124].includes(type)) throw new Error("model.wasm has an unsupported global type");\n          if (bytes[pos++] !== 0) throw new Error("model.wasm must not define mutable globals");\n          const op = bytes[pos++];\n          if (op === 65 || op === 66) {\n            let n = 0;\n            while (pos < end && n++ < 10 && bytes[pos++] & 128) {\n            }\n          } else if (op === 67) pos += 4;\n          else if (op === 68) pos += 8;\n          else throw new Error("model.wasm has an unsupported global initializer");\n          if (pos >= end || bytes[pos++] !== 11) throw new Error("model.wasm has an invalid global initializer");\n        }\n      }\n      pos = end;\n    }\n  }\n  function checkModel(modelModule, mathlibModule) {\n    const math = new Set(\n      WebAssembly.Module.exports(mathlibModule).filter((e) => e.kind === "function").map((e) => e.name)\n    );\n    let memory = false;\n    for (const { module, name, kind } of WebAssembly.Module.imports(modelModule)) {\n      if (module === "env" && name === "memory" && kind === "memory") memory = true;\n      else if (module === "math" && kind === "function" && math.has(name)) continue;\n      else {\n        throw new Error(\n          `model.wasm imports ${module}.${name}, which is neither env.memory nor a function of mathlib`\n        );\n      }\n    }\n    if (!memory) throw new Error("model.wasm does not import env.memory");\n  }\n  function jsonText(value) {\n    return typeof value === "string" ? value : JSON.stringify(value);\n  }\n  async function createRuntime({\n    engine,\n    mathlib,\n    console: out = globalThis.console,\n    maxMemoryPages = MAX_MEMORY_PAGES\n  }) {\n    const engineModule = await toModule(engine, "engine.wasm");\n    const mathlibModule = await toModule(mathlib, "mathlib.wasm");\n    checkRuntime(engineModule, mathlibModule);\n    const minimum = memoryMinimum(engineModule);\n    if (!Number.isSafeInteger(maxMemoryPages) || maxMemoryPages < minimum || maxMemoryPages > MAX_MEMORY_PAGES) {\n      throw new Error(`maxMemoryPages must be between ${minimum} and ${MAX_MEMORY_PAGES}`);\n    }\n    return {\n      /** The compiled modules, which a worker gets to create a runtime of its own. */\n      modules: { engine: engineModule, mathlib: mathlibModule },\n      /**\n       * Instantiates a model (a `WebAssembly.Module` or the bytes of `model.wasm`) with its\n       * description (the text of `model.json` or its object) and loads it into a new engine.\n       * Asynchronous (a Promise of the instance): bytes compile with `WebAssembly.compile`, which\n       * browsers require for modules above a few KB on the main thread.\n       */\n      async instantiate(modelWasm, description) {\n        const descText = jsonText(description);\n        const desc = typeof description === "string" ? JSON.parse(description) : description;\n        if (desc.contract !== MODEL_CONTRACT) {\n          throw new Error(`the description is of the contract ${desc.contract}, not ${MODEL_CONTRACT}`);\n        }\n        if (!(modelWasm instanceof WebAssembly.Module)) validateModelBytes(modelWasm);\n        const modelModule = await toModule(modelWasm, "model.wasm");\n        const history = WebAssembly.Module.customSections(modelModule, "sbmlsim.history");\n        const nDelays = desc.n_delays ?? 0;\n        if (history.length > 1 || history.length === 0 && nDelays !== 0) {\n          throw new Error("model.wasm history sites differ from its description");\n        }\n        if (history.length === 1) {\n          const info = new DataView(history[0]);\n          if (info.byteLength !== 8 || info.getUint32(0, true) !== 1 || info.getUint32(4, true) !== nDelays) {\n            throw new Error("model.wasm history metadata differs from its description");\n          }\n        }\n        checkModel(modelModule, mathlibModule);\n        const memory = new WebAssembly.Memory({ initial: minimum, maximum: maxMemoryPages });\n        const wasi = wasiImports(engineModule, () => memory, out);\n        const math = await WebAssembly.instantiate(mathlibModule, {});\n        const model = await WebAssembly.instantiate(modelModule, {\n          env: { memory },\n          math: math.exports\n        });\n        let engineInstance;\n        try {\n          engineInstance = await WebAssembly.instantiate(engineModule, {\n            env: { memory },\n            model: model.exports,\n            [WASI]: wasi.imports\n          });\n        } catch (error) {\n          if (error instanceof WebAssembly.LinkError) {\n            throw new Error(`engine.wasm does not instantiate with the model: ${error.message}`);\n          }\n          throw error;\n        }\n        const instance = new Instance(engineInstance.exports, memory, wasi.calls);\n        if (engineInstance.exports._initialize) {\n          instance.guard(() => engineInstance.exports._initialize());\n        }\n        instance.initializeCalls = wasi.calls();\n        instance.load(descText, desc.id);\n        return instance;\n      }\n    };\n  }\n  var preparedStates = /* @__PURE__ */ new WeakMap();\n  var PreparedSimulation = class {\n    constructor(token, state) {\n      if (token !== preparedStates) throw new TypeError("use instance.prepare()");\n      preparedStates.set(this, state);\n    }\n    get nValues() {\n      return preparedStates.get(this).nValues;\n    }\n    close() {\n      const state = preparedStates.get(this);\n      if (state.closed) return;\n      state.closed = true;\n      state.instance.plans.delete(this);\n      if (!state.instance.trapped) state.instance.guard(() => state.instance.exports.engine_release_prepared(state.handle));\n    }\n  };\n  var Instance = class {\n    constructor(exports, memory, calls) {\n      this.exports = exports;\n      this.memory = memory;\n      this.calls = calls;\n      this.trapped = false;\n      this.encoder = new TextEncoder();\n      this.decoder = new TextDecoder();\n      this.inputPtr = 0;\n      this.inputCapacity = 0;\n      this.plans = /* @__PURE__ */ new Set();\n    }\n    /** The number of calls of the engine into WASI since the instance was created. */\n    wasiCalls() {\n      return this.calls();\n    }\n    /** The size of the memory of the instance in bytes. */\n    memoryBytes() {\n      return this.memory.buffer.byteLength;\n    }\n    /**\n     * The counters of the integrator of the last simulation: an object of the names of\n     * `STATS_NAMES` to numbers (`engine_stats_ptr`).\n     */\n    stats() {\n      if (this.trapped) {\n        throw new EngineError("internal", "the instance trapped before and is not used again");\n      }\n      const ptr = this.guard(() => this.exports.engine_stats_ptr()) >>> 0;\n      const view = new DataView(this.memory.buffer, ptr, 8 * STATS_NAMES.length);\n      return Object.fromEntries(STATS_NAMES.map((name, i) => [name, Number(view.getBigUint64(8 * i, true))]));\n    }\n    /** The address of the output of the last call in the memory (aligned to 8). */\n    outputAddress() {\n      return this.exports.engine_output_ptr() >>> 0;\n    }\n    /**\n     * Simulates a request (an object of the contract `simulation/1` or its JSON text); gives\n     * `{ columns, nRows, values }`, the values a `Float64Array` of the rows, row major. The result\n     * does not depend on the simulations before. An error of the engine is an `EngineError`.\n     */\n    simulate(request) {\n      const code = this.call(this.exports.engine_simulate, jsonText(request));\n      const output = this.output();\n      if (code === 0 && output.table) return output.table;\n      if (code === 1 && output.error) throw output.error;\n      throw new EngineError("internal", `engine_simulate returned ${code} with the output ${JSON.stringify(output.header)}`);\n    }\n    /** Prepares fixed request structure once; varying IDs define the raw value order. */\n    prepare(request, { varying = [] } = {}) {\n      if (!["engine_prepare", "engine_simulate_prepared", "engine_release_prepared"].every((name) => typeof this.exports[name] === "function")) {\n        throw new UnsupportedError("this engine does not support prepared requests");\n      }\n      if (this.trapped) throw new EngineError("internal", "the instance trapped before and is not used again");\n      if (this.plans.size >= 32) throw new EngineError("request", "prepared request limit reached (32)");\n      const code = this.call(this.exports.engine_prepare, jsonText({ request: typeof request === "string" ? JSON.parse(request) : request, varying }));\n      const output = this.output();\n      if (code !== 0) throw output.error ?? new EngineError("internal", "engine_prepare failed without an error");\n      const { prepared: handle2, n_values: nValues } = output.header;\n      if (!Number.isInteger(handle2) || !Number.isSafeInteger(nValues) || nValues < 0) throw new EngineError("internal", "invalid prepared request output");\n      const plan = new PreparedSimulation(preparedStates, { instance: this, handle: handle2, nValues, closed: false });\n      this.plans.add(plan);\n      return plan;\n    }\n    simulatePrepared(plan, values = []) {\n      const state = preparedStates.get(plan);\n      if (!state || state.instance !== this || state.closed) throw new EngineError("request", "foreign or closed prepared request");\n      if (this.trapped) throw new EngineError("internal", "the instance trapped before and is not used again");\n      if (!Array.isArray(values) && !(values instanceof Float64Array) || values.length !== state.nValues || values.findIndex((value) => !Number.isFinite(value)) >= 0) {\n        throw new EngineError("request", `prepared request requires ${state.nValues} finite values`);\n      }\n      const ptr = this.inputBuffer(values.length * 8);\n      const view = new DataView(this.memory.buffer, ptr, values.length * 8);\n      for (let i = 0; i < values.length; i++) view.setFloat64(i * 8, values[i], true);\n      const code = this.guard(() => this.exports.engine_simulate_prepared(state.handle, ptr, values.length));\n      const output = this.output();\n      if (code === 0 && output.table) return output.table;\n      if (output.error) throw output.error;\n      throw new EngineError("internal", "engine_simulate_prepared returned invalid output");\n    }\n    inputBuffer(length) {\n      if (length > this.inputCapacity) {\n        const capacity = Math.max(length, this.inputCapacity * 2, 256);\n        const ptr = this.guard(() => this.exports.engine_alloc(capacity)) >>> 0;\n        if (!ptr) throw new EngineError("internal", "engine_alloc found no memory for the input");\n        if (this.inputCapacity) this.guard(() => this.exports.engine_free(this.inputPtr, this.inputCapacity));\n        this.inputPtr = ptr;\n        this.inputCapacity = capacity;\n      }\n      return this.inputPtr;\n    }\n    /** Loads the description (`engine_load`). */\n    load(descText, id) {\n      const before = this.calls();\n      const code = this.call(this.exports.engine_load, descText);\n      this.loadCalls = this.calls() - before;\n      if (code !== 0) {\n        const { error } = this.output();\n        const message = error?.kind === "load" ? error.detail : `engine_load failed: ${error?.message}`;\n        throw new Error(`the engine refuses the description of model ${id}: ${message}`);\n      }\n    }\n    /** Calls `engine_load` or `engine_simulate` with `text` in a buffer of the engine. */\n    call(fn, text) {\n      if (this.trapped) {\n        throw new EngineError("internal", "the instance trapped before and is not used again");\n      }\n      const bytes = this.encoder.encode(text);\n      const ptr = this.inputBuffer(bytes.length);\n      new Uint8Array(this.memory.buffer, ptr, bytes.length).set(bytes);\n      const code = this.guard(() => fn(ptr, bytes.length));\n      return code;\n    }\n    /**\n     * Runs a call of the engine; whatever it throws marks the instance and is an `EngineError`\n     * `internal` ([`trapError`]), as a trap of the wasmtime host, but `WasiExit`.\n     */\n    guard(call) {\n      try {\n        return call();\n      } catch (error) {\n        this.trapped = true;\n        if (error instanceof WasiExit) throw error;\n        let panic = null;\n        try {\n          const { error: e } = this.output();\n          if (e?.engineKind === "panic") panic = e.detail;\n        } catch {\n        }\n        throw trapError(error, panic);\n      }\n    }\n    /** The output of the last call: `{ header, table }` or `{ header, error }`. */\n    output() {\n      const ptr = this.exports.engine_output_ptr() >>> 0;\n      const len = this.exports.engine_output_len() >>> 0;\n      const buffer = this.memory.buffer;\n      if (len < 4 || ptr + len > buffer.byteLength) {\n        throw new EngineError("internal", "the output of the engine is outside of the memory");\n      }\n      const view = new DataView(buffer, ptr, len);\n      const headerLength = view.getUint32(0, true);\n      if (4 + headerLength > len) {\n        throw new EngineError("internal", "the header of the output of the engine is cut off");\n      }\n      const header = JSON.parse(this.decoder.decode(new Uint8Array(buffer, ptr + 4, headerLength)));\n      if (header.error) {\n        const { kind, message, flag, time } = header.error;\n        const known = ["load", "request", "unsupported", "integration"];\n        const error = known.includes(kind) ? new EngineError(kind, message, { flag, time }) : new EngineError("internal", kind === "panic" ? `engine.wasm panicked: ${message}` : `an error of the unknown kind ${kind}: ${message}`);\n        error.detail = message;\n        error.engineKind = kind;\n        return { header, error };\n      }\n      if (header.prepared !== void 0) return { header };\n      const { columns, n_rows: nRows } = header;\n      const start = 4 + headerLength + 7 & ~7;\n      const n = nRows * columns.length;\n      if (start + 8 * n !== len) {\n        throw new EngineError(\n          "internal",\n          `the output of the engine does not hold ${nRows} rows of ${columns.length} columns`\n        );\n      }\n      let values;\n      if (LITTLE_ENDIAN) {\n        values = new Float64Array(buffer, ptr + start, n).slice();\n      } else {\n        values = new Float64Array(n);\n        for (let i = 0; i < n; i += 1) values[i] = view.getFloat64(start + 8 * i, true);\n      }\n      return { header, table: { columns, nRows, values } };\n    }\n  };\n  function trapError(error, panic = null) {\n    const what = error instanceof WebAssembly.RuntimeError ? "trapped" : `failed (${error?.name ?? "Error"})`;\n    const message = `engine.wasm ${what}: ${error?.message ?? String(error)}${panic ? ` (${panic})` : ""}`;\n    const mapped = new EngineError("internal", message);\n    mapped.cause = error;\n    return mapped;\n  }\n\n  // worker.mjs\n  var runtime = null;\n  var models = /* @__PURE__ */ new Map();\n  var post = null;\n  var queue = Promise.resolve();\n  function errorFields(error) {\n    const fields = { name: error?.name ?? "Error", message: error?.message ?? String(error) };\n    for (const name of ["kind", "detail", "flag", "time"]) {\n      if (error?.[name] !== void 0) fields[name] = error[name];\n    }\n    return fields;\n  }\n  async function instanceOf(key) {\n    const model = models.get(key);\n    if (!model) throw new Error(`the pool posted no model ${key} to this worker`);\n    if (!runtime) throw new Error("the pool posted no runtime to this worker");\n    if (!model.instance) model.instance = (await runtime).instantiate(model.module, model.description);\n    try {\n      return await model.instance;\n    } catch (error) {\n      model.instance = null;\n      throw error;\n    }\n  }\n  async function handle(message) {\n    switch (message.type) {\n      case "init":\n        runtime = createRuntime({ engine: message.engine, mathlib: message.mathlib, maxMemoryPages: message.maxMemoryPages });\n        runtime.catch(() => {\n        });\n        return;\n      case "model":\n        models.set(message.key, { module: message.module, description: message.description, instance: null, plans: /* @__PURE__ */ new Map(), definitions: /* @__PURE__ */ new Map() });\n        return;\n      case "releasePrepared": {\n        const model = models.get(message.key);\n        model?.plans.get(message.planId)?.close();\n        model?.plans.delete(message.planId);\n        model?.definitions.delete(message.planId);\n        return;\n      }\n      case "removeModel":\n        for (const plan of models.get(message.key)?.plans.values() ?? []) plan.close();\n        models.delete(message.key);\n        return;\n      case "prepare":\n      case "simulatePrepared":\n      case "simulate": {\n        let instance = null;\n        try {\n          instance = await instanceOf(message.key);\n          const model = models.get(message.key);\n          let plan;\n          if (message.type !== "simulate") {\n            if (message.request !== void 0) model.definitions.set(message.planId, { request: message.request, varying: message.varying });\n            plan = model.plans.get(message.planId);\n            if (!plan) {\n              const definition = model.definitions.get(message.planId);\n              if (!definition) throw new Error("missing prepared request definition");\n              plan = instance.prepare(definition.request, { varying: definition.varying });\n              model.plans.set(message.planId, plan);\n            }\n          }\n          if (message.type === "prepare") {\n            post({ type: "prepared", id: message.id, nValues: plan.nValues });\n            return;\n          }\n          const { columns, nRows, values } = plan ? instance.simulatePrepared(plan, message.values) : instance.simulate(message.request);\n          post({ type: "result", id: message.id, columns, nRows, values }, [values.buffer]);\n        } catch (error) {\n          if (instance?.trapped) {\n            const model = models.get(message.key);\n            model.instance = null;\n            model.plans.clear();\n          }\n          post({ type: "error", id: message.id, error: errorFields(error) });\n        }\n        return;\n      }\n      default:\n        throw new Error(`a message of the unknown type ${message.type}`);\n    }\n  }\n  function receive(message) {\n    queue = queue.then(() => handle(message)).catch((error) => console.error("sbmlsim-core worker:", error));\n  }\n  if (typeof WorkerGlobalScope !== "undefined" && globalThis instanceof WorkerGlobalScope) {\n    post = (message, transfer = []) => globalThis.postMessage(message, transfer);\n    globalThis.addEventListener("message", (event) => receive(event.data));\n  } else {\n    const specifier = "node:worker_threads";\n    import(specifier).then(({ parentPort }) => {\n      if (!parentPort) throw new Error("worker.mjs runs as a Web Worker or as a worker thread of Node");\n      post = (message, transfer = []) => parentPort.postMessage(message, transfer);\n      parentPort.on("message", receive);\n    });\n  }\n})();\n';
  var ORIGIN = "https://sbmlsim-core.invalid";
  var BUNDLE_URL = `${ORIGIN}/bundle/bundle.json`;
  function blobWorkers(source = WORKER_SOURCE) {
    let url = null;
    return () => {
      if (typeof globalThis.Worker !== "function") throw new Error("this host has no Web Workers");
      url ??= URL.createObjectURL(new Blob([source], { type: "text/javascript" }));
      return new globalThis.Worker(url);
    };
  }
  function respond(body, status = 200) {
    return new Response(body, { status });
  }
  function createInlineEndpoint({ bundle, files = {}, engine, mathlib, workers, createWorker = blobWorkers() }) {
    if (bundle === null || typeof bundle !== "object") throw new Error("createInlineEndpoint needs the bundle as an object");
    if (!engine || !mathlib) throw new Error("createInlineEndpoint needs the engine and the math library");
    const text = JSON.stringify(bundle);
    const base = new URL(BUNDLE_URL);
    const load = async (url) => {
      const target = new URL(String(url), base);
      if (target.href === base.href) return respond(text);
      const prefix = new URL("./", base).href;
      if (!target.href.startsWith(prefix)) return respond(`${target.href} is no file of the bundle`, 404);
      const file = files[decodeURIComponent(target.href.slice(prefix.length))];
      if (file === void 0) return respond(`${target.href} is no file of the bundle`, 404);
      return respond(file);
    };
    return createPanelEndpoint({
      bundle: BUNDLE_URL,
      runtime: { engine, mathlib },
      ...workers === void 0 ? {} : { workers },
      createWorker,
      fetch: load
    });
  }
  return __toCommonJS(inline_exports);
})();
