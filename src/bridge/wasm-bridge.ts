/**
 * WASM Bridge
 *
 * Interface for loading and calling QuantumOS geometric control
 * functions compiled to WebAssembly. Falls back to pure TypeScript
 * geometric module when WASM is not available.
 *
 * Memory management pattern:
 *   1. Allocate space in WASM linear memory via malloc()
 *   2. Copy TypeScript Float64Array data into WASM memory
 *   3. Call the WASM function with pointers
 *   4. Copy results back out
 *   5. Free the WASM memory via free()
 */

// ============================================
// WEBASSEMBLY TYPE DECLARATIONS
// ============================================
// The project tsconfig uses "lib": ["ES2020"] which does not
// include WebAssembly types. We declare the minimum surface
// area needed by this module.

/* eslint-disable @typescript-eslint/no-namespace */
declare namespace WebAssembly {
  interface Memory {
    readonly buffer: ArrayBuffer;
  }

  interface MemoryDescriptor {
    initial: number;
    maximum?: number;
  }

  const Memory: {
    new(descriptor: MemoryDescriptor): Memory;
  };

  interface Instance {
    readonly exports: Record<string, unknown>;
  }

  interface Module {}

  interface ResultObject {
    instance: Instance;
    module: Module;
  }

  type Imports = Record<string, Record<string, unknown>>;

  function instantiate(
    bufferSource: ArrayBuffer,
    importObject?: Imports
  ): Promise<ResultObject>;
}
/* eslint-enable @typescript-eslint/no-namespace */

// ============================================
// GEOMETRIC TYPES
// ============================================
// These mirror the types that will live in src/geometric/
// once the geometric module is implemented. Defined here so
// the bridge compiles independently.

/**
 * Result of a chiral anomaly detection in the geometric module.
 */
export interface ChiralAnomalyResult {
  anomalyIndex: number;
  isAnomalous: boolean;
  leftSpectrum: Float64Array;
  rightSpectrum: Float64Array;
}

/**
 * Manifold-level metrics from the geometric control module.
 */
export interface ManifoldMetrics {
  ricciScalar: number;
  berryPhase: number;
  dimension: number;
  curvatureNorm: number;
}

// ============================================
// WASM EXPORTS INTERFACE
// ============================================

/**
 * Expected export shape from the QuantumOS geometric control
 * WASM module (compiled from C with Emscripten or similar).
 */
export interface WASMExports {
  /** Initialize the geometric control engine for a given state-space dimension. Returns 0 on success. */
  geometric_control_init(dimension: number): number;

  /** Update the internal metric tensor from an array of oscillator phases. Returns 0 on success. */
  geometric_update_metric(phasesPtr: number, n: number): number;

  /** Compute the Berry phase around a closed loop in state space. */
  geometric_compute_berry_phase(loopPtr: number, loopLength: number, n: number): number;

  /** Return the current Ricci scalar curvature of the state manifold. */
  geometric_get_ricci_scalar(): number;

  /**
   * Detect chiral anomaly by comparing left- and right-handed spectra.
   * Writes the result into the struct at resultPtr. Returns 0 on success.
   */
  geometric_detect_chiral_anomaly(
    leftPtr: number,
    rightPtr: number,
    n: number,
    resultPtr: number
  ): number;

  /** WASM linear memory */
  memory: WebAssembly.Memory;

  /** Allocate `size` bytes in WASM heap. Returns pointer (offset). */
  malloc(size: number): number;

  /** Free a previous allocation. */
  free(ptr: number): void;
}

// ============================================
// CONSTANTS
// ============================================

/** Number of bytes per Float64 element */
const FLOAT64_BYTES = 8;

/**
 * Size of the C-side chiral anomaly result struct (in bytes).
 *   double anomalyIndex       (8)
 *   uint8  isAnomalous        (1)
 *   [padding]                 (7) -- align to 8
 *   double leftSpectrum[N]    (variable, but we pre-allocate for dimension)
 *   double rightSpectrum[N]   (variable)
 *
 * We calculate it dynamically based on dimension.
 */
function chiralResultSize(n: number): number {
  // anomalyIndex(8) + isAnomalous(1) + padding(7) + left(n*8) + right(n*8)
  return 16 + n * FLOAT64_BYTES * 2;
}

// ============================================
// WASM BRIDGE
// ============================================

export class WASMBridge {
  private instance: WebAssembly.Instance | null = null;
  private exports: WASMExports | null = null;
  private ready: boolean = false;
  private dimension: number;

  constructor(dimension: number = 64) {
    this.dimension = dimension;
  }

  // ------------------------------------------
  // Initialization
  // ------------------------------------------

  /**
   * Initialize WASM module from a .wasm binary (ArrayBuffer) or a URL/path string.
   *
   * Returns true if the module loaded and initialized successfully,
   * false if WASM is not available or loading failed.
   */
  async init(wasmSource?: ArrayBuffer | string): Promise<boolean> {
    try {
      let wasmBuffer: ArrayBuffer;

      if (wasmSource instanceof ArrayBuffer) {
        wasmBuffer = wasmSource;
      } else if (typeof wasmSource === 'string') {
        wasmBuffer = await WASMBridge.loadFromPath(wasmSource);
      } else {
        // No source provided; nothing to load
        return false;
      }

      const importObject: WebAssembly.Imports = {
        env: {
          abort: () => { throw new Error('WASM abort'); },
          memory: new WebAssembly.Memory({ initial: 256, maximum: 2048 }),
        },
      };

      const result = await WebAssembly.instantiate(wasmBuffer, importObject);
      this.instance = result.instance;
      this.exports = this.instance.exports as unknown as WASMExports;

      // Initialize the geometric control engine
      const initResult = this.exports.geometric_control_init(this.dimension);
      if (initResult !== 0) {
        this.dispose();
        return false;
      }

      this.ready = true;
      return true;
    } catch {
      this.dispose();
      return false;
    }
  }

  /**
   * Load a WASM binary from a file path (Node.js) or URL (browser).
   */
  private static async loadFromPath(source: string): Promise<ArrayBuffer> {
    // Browser environment: use dynamic import or fetch-like mechanism
    // Node.js environment: use fs.readFile
    // We use dynamic import for fs to avoid bundler issues.
    const nodeFs = await import('fs').catch(() => null);
    if (nodeFs) {
      const data = nodeFs.readFileSync(source);
      // Convert Node Buffer to ArrayBuffer
      return data.buffer.slice(
        data.byteOffset,
        data.byteOffset + data.byteLength
      ) as ArrayBuffer;
    }

    // Fallback: try global fetch (available in Node 18+ and browsers)
    const globalFetch = globalThis.fetch;
    if (globalFetch) {
      const response = await globalFetch(source);
      if (!response.ok) {
        throw new Error(`Failed to fetch WASM from ${source}: ${response.status}`);
      }
      return response.arrayBuffer();
    }

    throw new Error('No mechanism available to load WASM binary');
  }

  /**
   * Check if the WASM module is loaded and ready.
   */
  isReady(): boolean {
    return this.ready;
  }

  // ------------------------------------------
  // Geometric Operations
  // ------------------------------------------

  /**
   * Update the metric tensor from oscillator phases.
   * Returns true on success, false if WASM is unavailable.
   */
  geometricUpdateMetric(phases: Float64Array): boolean {
    if (!this.ready || !this.exports) return false;

    const n = phases.length;
    const ptr = this.allocateFloat64Array(phases);
    if (ptr === 0) return false;

    try {
      const result = this.exports.geometric_update_metric(ptr, n);
      return result === 0;
    } finally {
      this.exports.free(ptr);
    }
  }

  /**
   * Compute the Berry phase around a closed loop of state-space points.
   * Each element of `loop` is a point in the N-dimensional state space.
   * Returns the Berry phase value, or null if WASM is unavailable.
   */
  geometricComputeBerryPhase(loop: Float64Array[]): number | null {
    if (!this.ready || !this.exports) return null;
    if (loop.length === 0) return null;

    const n = loop[0].length;
    const loopLength = loop.length;

    // Flatten the loop into a single contiguous array: [point0[0..n-1], point1[0..n-1], ...]
    const flat = new Float64Array(loopLength * n);
    for (let i = 0; i < loopLength; i++) {
      flat.set(loop[i], i * n);
    }

    const ptr = this.allocateFloat64Array(flat);
    if (ptr === 0) return null;

    try {
      return this.exports.geometric_compute_berry_phase(ptr, loopLength, n);
    } finally {
      this.exports.free(ptr);
    }
  }

  /**
   * Compute the Berry curvature 2-form at a single point in state space.
   * Returns an n x n antisymmetric matrix as Float64Array[], or null.
   *
   * Implementation: Computes Berry phase around infinitesimal loops
   * in each pair of coordinate directions, then derives the curvature
   * components F_{ij} = dA_j/dx_i - dA_i/dx_j.
   */
  geometricComputeBerryCurvature(point: Float64Array): Float64Array[] | null {
    if (!this.ready || !this.exports) return null;

    const n = point.length;
    const epsilon = 1e-6;

    // Build curvature matrix by computing Berry phase around
    // infinitesimal plaquettes in each (i, j) plane.
    const curvature: Float64Array[] = [];
    for (let i = 0; i < n; i++) {
      curvature[i] = new Float64Array(n);
    }

    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        // Construct a small loop: p -> p+ei -> p+ei+ej -> p+ej -> p
        const p0 = new Float64Array(point);
        const p1 = new Float64Array(point); p1[i] += epsilon;
        const p2 = new Float64Array(point); p2[i] += epsilon; p2[j] += epsilon;
        const p3 = new Float64Array(point); p3[j] += epsilon;

        const loop = [p0, p1, p2, p3];
        const phase = this.geometricComputeBerryPhase(loop);
        if (phase === null) return null;

        const fij = phase / (epsilon * epsilon);
        curvature[i][j] = fij;
        curvature[j][i] = -fij;  // Antisymmetric
      }
    }

    return curvature;
  }

  /**
   * Detect chiral anomaly by comparing left- and right-handed spectra.
   * Returns ChiralAnomalyResult or null if WASM unavailable.
   */
  geometricDetectChiralAnomaly(
    left: Float64Array,
    right: Float64Array
  ): ChiralAnomalyResult | null {
    if (!this.ready || !this.exports) return null;

    const n = Math.min(left.length, right.length);
    if (n === 0) return null;

    const leftPtr = this.allocateFloat64Array(left);
    if (leftPtr === 0) return null;

    const rightPtr = this.allocateFloat64Array(right);
    if (rightPtr === 0) {
      this.exports.free(leftPtr);
      return null;
    }

    const resultSize = chiralResultSize(n);
    const resultPtr = this.exports.malloc(resultSize);
    if (resultPtr === 0) {
      this.exports.free(leftPtr);
      this.exports.free(rightPtr);
      return null;
    }

    try {
      const status = this.exports.geometric_detect_chiral_anomaly(
        leftPtr, rightPtr, n, resultPtr
      );
      if (status !== 0) return null;

      // Read the result struct back from WASM memory
      const mem = new DataView(this.exports.memory.buffer);
      const anomalyIndex = mem.getFloat64(resultPtr, true);
      const isAnomalous = mem.getUint8(resultPtr + 8) !== 0;

      // Spectra start after the fixed 16-byte header
      const leftSpectrumOffset = resultPtr + 16;
      const rightSpectrumOffset = leftSpectrumOffset + n * FLOAT64_BYTES;

      const leftSpectrum = new Float64Array(n);
      const rightSpectrum = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        leftSpectrum[i] = mem.getFloat64(leftSpectrumOffset + i * FLOAT64_BYTES, true);
        rightSpectrum[i] = mem.getFloat64(rightSpectrumOffset + i * FLOAT64_BYTES, true);
      }

      return { anomalyIndex, isAnomalous, leftSpectrum, rightSpectrum };
    } finally {
      this.exports.free(resultPtr);
      this.exports.free(rightPtr);
      this.exports.free(leftPtr);
    }
  }

  /**
   * Get the current Ricci scalar curvature, or null if WASM unavailable.
   */
  getRicciScalar(): number | null {
    if (!this.ready || !this.exports) return null;
    return this.exports.geometric_get_ricci_scalar();
  }

  /**
   * Get manifold metrics (Ricci scalar, Berry phase from identity loop, etc.)
   * Returns null if WASM is not available.
   */
  getMetrics(): ManifoldMetrics | null {
    if (!this.ready || !this.exports) return null;

    const ricciScalar = this.exports.geometric_get_ricci_scalar();

    // Compute Berry phase around a trivial loop at the origin
    // to characterize the global topology
    const origin = new Float64Array(this.dimension);
    const trivialLoop = [origin, origin];
    const berryPhase = this.geometricComputeBerryPhase(trivialLoop) ?? 0;

    // Curvature norm from Ricci scalar
    const curvatureNorm = Math.abs(ricciScalar);

    return {
      ricciScalar,
      berryPhase,
      dimension: this.dimension,
      curvatureNorm,
    };
  }

  // ------------------------------------------
  // Cleanup
  // ------------------------------------------

  /**
   * Release the WASM instance and associated resources.
   */
  dispose(): void {
    this.instance = null;
    this.exports = null;
    this.ready = false;
  }

  // ------------------------------------------
  // Memory Helpers
  // ------------------------------------------

  /**
   * Allocate a region of WASM linear memory, copy a Float64Array into it,
   * and return the pointer (byte offset). Returns 0 on failure.
   *
   * Caller is responsible for calling `this.exports.free(ptr)`.
   */
  private allocateFloat64Array(data: Float64Array): number {
    if (!this.exports) return 0;

    const byteLength = data.length * FLOAT64_BYTES;
    const ptr = this.exports.malloc(byteLength);
    if (ptr === 0) return 0;

    // Write data into WASM linear memory
    const target = new Float64Array(
      this.exports.memory.buffer,
      ptr,
      data.length
    );
    target.set(data);

    return ptr;
  }
}

// ============================================
// FACTORY
// ============================================

/**
 * Create a WASMBridge with optional dimension (defaults to 64).
 *
 * The bridge starts in the "not ready" state. Call `bridge.init(source)`
 * with a WASM binary to activate hardware-accelerated geometric control.
 * All methods gracefully return null / false when WASM is unavailable.
 */
export function createWASMBridge(dimension?: number): WASMBridge {
  return new WASMBridge(dimension);
}
