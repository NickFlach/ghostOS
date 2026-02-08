/**
 * ghostOS Resonance Engine
 *
 * The core dynamics module implementing the fundamental equation:
 *   dx/dt = f(x) - λx
 *
 * Where resonance emerges in the balanced regime between
 * divergence (λ→0) and collapse (λ→∞).
 */

import {
  DEFAULT_CONFIG,
  DAMPING_DEFAULT,
  DAMPING_MIN,
  DAMPING_MAX,
  RESONANCE_BAND,
  PHI,
  EMERGENCE_COEFFICIENT,
  TAU,
  CHIRAL_ETA_DEFAULT,
} from '../constants';
import type { Signal } from '../signal';
import { ChiralEngine, type ChiralState, type ChiralConfig } from '../chiral';
import { OscillatorManifold, type ManifoldMetrics } from '../geometric/manifold';
import { createBerryTracker, type BerryState } from '../geometric/berry';
import { detectChiralAnomaly, type ChiralAnomalyResult } from '../geometric/anomaly';

// ============================================
// TYPES
// ============================================

export interface ResonanceConfig {
  lambda?: number;
  adaptiveRate?: number;
  targetCoherence?: number;
  dt?: number;
  /** Enable chiral stability mechanisms */
  enableChiral?: boolean;
  /** Chiral configuration options */
  chiralConfig?: ChiralConfig;
}

export interface ResonantState {
  state: Float64Array;
  coherence: number;
  energy: number;
  dominantFreq: number;
  isStable: boolean;
  emergence: Float64Array;
  phase: number;
  lambda: number;
  /** Chiral state (if chiral mode enabled) */
  chiral?: ChiralState;
  /** Geometric manifold metrics (curvature, Berry phase) */
  geometric?: {
    manifold: ManifoldMetrics;
    berry: BerryState;
    anomaly: ChiralAnomalyResult | null;
  };
}

export interface ModeInfo {
  frequency: number;
  amplitude: number;
  phase: number;
  stability: number;
}

// ============================================
// RESONANCE ENGINE
// ============================================

export class ResonanceEngine {
  private lambda: number;
  private adaptiveRate: number;
  private targetCoherence: number;
  private dt: number;
  private state: Float64Array;
  private previousState: Float64Array;
  private velocity: Float64Array;  // For chiral wave dynamics
  private phase: number = 0;
  private tick: number = 0;

  // Chiral stability components
  private enableChiral: boolean;
  private chiralEngine: ChiralEngine | null = null;

  // Geometric control components
  private manifold: OscillatorManifold;
  private berryTracker: ReturnType<typeof createBerryTracker>;
  private lastAnomaly: ChiralAnomalyResult | null = null;

  constructor(config: ResonanceConfig = {}) {
    this.lambda = config.lambda ?? DAMPING_DEFAULT;
    this.adaptiveRate = config.adaptiveRate ?? DEFAULT_CONFIG.adaptiveRate;
    this.targetCoherence = config.targetCoherence ?? DEFAULT_CONFIG.targetCoherence;
    this.dt = config.dt ?? DEFAULT_CONFIG.dt;
    this.state = new Float64Array(DEFAULT_CONFIG.dimensions);
    this.previousState = new Float64Array(DEFAULT_CONFIG.dimensions);
    this.velocity = new Float64Array(DEFAULT_CONFIG.dimensions);

    // Initialize geometric control layer
    this.manifold = new OscillatorManifold(DEFAULT_CONFIG.dimensions);
    this.berryTracker = createBerryTracker();

    // Initialize chiral stability if enabled
    this.enableChiral = config.enableChiral ?? true;
    if (this.enableChiral) {
      this.chiralEngine = new ChiralEngine({
        eta: CHIRAL_ETA_DEFAULT,
        gamma: this.lambda,
        ...config.chiralConfig,
      });
    }
  }

  /**
   * Process a signal through the resonance layer
   * Applies: dx/dt = f(x) - λx
   *
   * With chiral enhancement:
   * - Non-reciprocal coupling for directional stability
   * - CISS coherence boost
   * - Topological protection of edge modes
   */
  process(signal: Signal): ResonantState {
    const n = signal.dimensions;

    // Ensure state arrays match signal dimensions
    if (this.state.length !== n) {
      this.state = new Float64Array(n);
      this.previousState = new Float64Array(n);
      this.velocity = new Float64Array(n);
    }

    // Store previous state for emergence calculation
    this.previousState.set(this.state);

    // Apply resonant dynamics with optional chiral enhancement
    if (this.enableChiral && this.chiralEngine) {
      this.processWithChiral(signal, n);
    } else {
      this.processStandard(signal, n);
    }

    // Calculate emergence (the residue of transformation)
    const emergence = this.calculateEmergence();

    // Calculate metrics
    const energy = this.calculateEnergy();
    let coherence = this.calculateCoherence();
    const dominantFreq = this.estimateDominantFrequency();

    // Apply CISS coherence enhancement if chiral is enabled
    if (this.enableChiral && this.chiralEngine) {
      coherence = this.chiralEngine.applyCISSEnhancement(coherence);
    }

    const isStable = this.checkStability(energy, coherence);

    // Update geometric control layer
    const gradient = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      gradient[i] = this.state[i] - this.previousState[i];
    }
    this.manifold.updateMetric(
      { coordinates: this.state, tangentVector: gradient, dimension: n },
      gradient
    );
    this.berryTracker.update(this.state, gradient);

    // Chiral anomaly detection from mode analysis
    if (this.enableChiral) {
      const modes = this.detectModes();
      if (modes.length >= 2) {
        const leftModes = new Float64Array(modes.filter((_, i) => i % 2 === 0).map(m => m.amplitude));
        const rightModes = new Float64Array(modes.filter((_, i) => i % 2 === 1).map(m => m.amplitude));
        this.lastAnomaly = detectChiralAnomaly(leftModes, rightModes);
      }
    }

    // Adaptive damping
    this.adaptDamping(coherence);

    // Update phase
    this.phase = (this.phase + dominantFreq * this.dt) % TAU;
    this.tick++;

    // Build result with optional chiral state
    const result: ResonantState = {
      state: new Float64Array(this.state),
      coherence,
      energy,
      dominantFreq,
      isStable,
      emergence,
      phase: this.phase,
      lambda: this.lambda,
    };

    if (this.enableChiral && this.chiralEngine) {
      result.chiral = this.chiralEngine.getState();
    }

    // Attach geometric state
    result.geometric = {
      manifold: this.manifold.getMetrics(),
      berry: this.berryTracker.getState(),
      anomaly: this.lastAnomaly,
    };

    return result;
  }

  /**
   * Standard resonance processing (original algorithm)
   */
  private processStandard(signal: Signal, n: number): void {
    for (let i = 0; i < n; i++) {
      const fx = signal.data[i] + this.nonlinearCoupling(i);
      const constraint = this.lambda * this.state[i];
      const derivative = fx - constraint;
      this.state[i] += derivative * this.dt;
    }
  }

  /**
   * Chiral-enhanced resonance processing
   * Implements non-reciprocal coupling and wave dynamics
   */
  private processWithChiral(signal: Signal, n: number): void {
    // First, add input signal to state
    for (let i = 0; i < n; i++) {
      this.state[i] += signal.data[i] * this.dt;
    }

    // Apply chiral dynamics (non-reciprocal sine-Gordon)
    const result = this.chiralEngine!.applyChiralDynamics(
      this.state,
      this.velocity,
      this.dt
    );

    this.state = result.state;
    this.velocity = result.velocity;

    // Apply chiral-aware nonlinear coupling
    for (let i = 0; i < n; i++) {
      const coupling = this.chiralNonlinearCoupling(i);
      const constraint = this.lambda * this.state[i];
      this.state[i] += (coupling - constraint) * this.dt;
    }

    // Apply chiral gain to stabilize preferred handedness
    const gainedState = this.chiralEngine!.applyChiralGain(this.state, 0.05);
    this.state = gainedState;
  }

  /**
   * Chiral-aware nonlinear coupling
   * Uses non-reciprocal coupling strengths J_{i,j} ≠ J_{j,i}
   */
  private chiralNonlinearCoupling(index: number): number {
    const n = this.state.length;
    let coupling = 0;

    for (let j = 0; j < n; j++) {
      if (j !== index) {
        // Get non-reciprocal coupling from chiral engine
        const chiralCoupling = this.chiralEngine!.calculateCoupling(
          index, j, PHI / n
        );

        // Use forward or backward coupling based on direction
        const strength = j > index ? chiralCoupling.forward : chiralCoupling.backward;
        coupling += strength * Math.sin(this.state[j] - this.state[index]);
      }
    }

    return coupling;
  }

  /**
   * Nonlinear coupling between state dimensions
   * Implements soft coupling inspired by Kuramoto oscillators
   */
  private nonlinearCoupling(index: number): number {
    const n = this.state.length;
    let coupling = 0;

    // Sum of sin(θj - θi) style coupling
    for (let j = 0; j < n; j++) {
      if (j !== index) {
        coupling += Math.sin(this.state[j] - this.state[index]);
      }
    }

    return (coupling * PHI) / n;  // Scale by golden ratio
  }

  /**
   * Calculate emergence residue (Ξ)
   * This is the difference between two transformation orderings
   */
  private calculateEmergence(): Float64Array {
    const n = this.state.length;
    const emergence = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      // Emergence = RG(x) - GR(x) approximation
      // R: rotation (use phase shift)
      // G: golden scaling

      const rotated = this.state[(i + 1) % n];
      const scaled = this.state[i] * PHI;

      const rg = rotated * PHI;  // Rotate then scale
      const gr = scaled;          // Scale then rotate (identity for diagonal)

      emergence[i] = EMERGENCE_COEFFICIENT * (rg - gr);
    }

    return emergence;
  }

  /**
   * Calculate system energy (L2 norm squared)
   */
  private calculateEnergy(): number {
    let energy = 0;
    for (let i = 0; i < this.state.length; i++) {
      energy += this.state[i] * this.state[i];
    }
    return energy;
  }

  /**
   * Calculate coherence (order parameter R)
   * Measures how synchronized the state dimensions are
   */
  private calculateCoherence(): number {
    const n = this.state.length;
    let sumCos = 0;
    let sumSin = 0;

    for (let i = 0; i < n; i++) {
      sumCos += Math.cos(this.state[i]);
      sumSin += Math.sin(this.state[i]);
    }

    const avgCos = sumCos / n;
    const avgSin = sumSin / n;

    return Math.sqrt(avgCos * avgCos + avgSin * avgSin);
  }

  /**
   * Estimate dominant frequency from state evolution
   */
  private estimateDominantFrequency(): number {
    let maxDiff = 0;

    for (let i = 0; i < this.state.length; i++) {
      const diff = Math.abs(this.state[i] - this.previousState[i]);
      if (diff > maxDiff) maxDiff = diff;
    }

    return maxDiff / this.dt;
  }

  /**
   * Check if system is in stable resonance
   */
  private checkStability(energy: number, coherence: number): boolean {
    // Stable if coherence in resonance band and energy bounded
    return (
      coherence >= RESONANCE_BAND.min &&
      coherence <= RESONANCE_BAND.max &&
      energy < 100 &&
      energy > 0.01
    );
  }

  /**
   * Adaptive damping control
   * Adjusts λ to maintain target coherence
   */
  private adaptDamping(currentCoherence: number): void {
    const error = this.targetCoherence - currentCoherence;

    // If coherence too low, reduce damping (allow more dynamics)
    // If coherence too high, increase damping (prevent lock-in)
    this.lambda -= error * this.adaptiveRate;

    // Clamp to safe range
    this.lambda = Math.max(DAMPING_MIN, Math.min(DAMPING_MAX, this.lambda));
  }

  /**
   * Manually adjust damping coefficient
   */
  adjustLambda(factor: number): void {
    this.lambda *= factor;
    this.lambda = Math.max(DAMPING_MIN, Math.min(DAMPING_MAX, this.lambda));
  }

  /**
   * Inject noise to break degeneracy
   */
  perturbState(intensity: number): void {
    for (let i = 0; i < this.state.length; i++) {
      this.state[i] += (Math.random() - 0.5) * 2 * intensity;
    }
  }

  /**
   * Reset state to initial conditions
   */
  reset(): void {
    this.state.fill(0);
    this.previousState.fill(0);
    this.velocity.fill(0);
    this.phase = 0;
    this.tick = 0;
    this.lambda = DAMPING_DEFAULT;
    if (this.chiralEngine) {
      this.chiralEngine.reset();
    }
    this.manifold.reset();
    this.berryTracker.reset();
    this.lastAnomaly = null;
  }

  /**
   * Get chiral stability status
   */
  getChiralState(): ChiralState | null {
    return this.chiralEngine?.getState() ?? null;
  }

  /**
   * Check if system is chirally stable
   */
  isChirallyStable(): boolean {
    if (!this.chiralEngine) return true;
    return this.chiralEngine.getStabilityRegime() === 'stable';
  }

  /**
   * Optimize chiral parameters for maximum stability
   */
  optimizeChiralStability(): void {
    if (this.chiralEngine) {
      this.chiralEngine.optimizeForStability();
    }
  }

  /**
   * Get current damping coefficient
   */
  getLambda(): number {
    return this.lambda;
  }

  /**
   * Detect resonant modes in the current state
   */
  detectModes(): ModeInfo[] {
    const modes: ModeInfo[] = [];
    const n = this.state.length;

    // Simple peak detection in frequency domain (approximation)
    for (let k = 1; k < n / 2; k++) {
      let realPart = 0;
      let imagPart = 0;

      for (let i = 0; i < n; i++) {
        const angle = (TAU * k * i) / n;
        realPart += this.state[i] * Math.cos(angle);
        imagPart += this.state[i] * Math.sin(angle);
      }

      const amplitude = Math.sqrt(realPart * realPart + imagPart * imagPart) / n;
      const phase = Math.atan2(imagPart, realPart);

      if (amplitude > 0.1) {
        modes.push({
          frequency: k,
          amplitude,
          phase,
          stability: amplitude > 0.3 ? 1 : amplitude / 0.3,
        });
      }
    }

    return modes.sort((a, b) => b.amplitude - a.amplitude);
  }
}

// ============================================
// CONVENIENCE EXPORTS
// ============================================

export const createEngine = (config?: ResonanceConfig) =>
  new ResonanceEngine(config);
