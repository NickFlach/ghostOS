/**
 * ghostmagic Chiral Dynamics Module
 *
 * Implements chiral stability mechanisms based on recent quantum research:
 *
 * 1. Non-reciprocal coupling (J_{i,j} ≠ J_{j,i}) for directional stability
 * 2. Time-reversal symmetry breaking through chiral gain
 * 3. Topological protection via chiral edge modes
 * 4. CISS (Chiral-Induced Spin Selectivity) for coherence enhancement
 *
 * Core equation (dissipative non-reciprocal sine-Gordon):
 *   φₜₜ - φₓₓ + sin(φ) = -ηφₓ - Γφₜ
 *
 * Where:
 *   η = non-reciprocity strength (breaks spatial inversion symmetry)
 *   Γ = on-site damping
 *   c = η/Γ determines propagation direction and stability
 *
 * Key insight: Waves propagating with c = η/Γ are stable,
 * while mirror-symmetric counterparts (c = -η/Γ) are unstable.
 * Chirality acts as a "master switch" for stability selection.
 *
 * References:
 * - arXiv:2505.02718 - Chiral Gain-Induced Time-Reversal Symmetry Breaking
 * - arXiv:2512.08092 - Nonreciprocity-Enabled Chiral Stability
 * - Nature Photonics 2025 - Chirality-induced quantum non-reciprocity
 */

import {
  PHI,
  PHI_INVERSE,
  TAU,
  CHIRAL_ETA_DEFAULT,
  CHIRAL_VELOCITY_THRESHOLD,
  NON_RECIPROCAL_ASYMMETRY,
  CHIRAL_STABILITY,
  TOPOLOGICAL_PROTECTION_FACTOR,
  CISS_COUPLING,
  NON_HERMITIAN,
  type Handedness,
} from '../constants';

// ============================================
// TYPES
// ============================================

export interface ChiralConfig {
  /** Non-reciprocity strength η */
  eta?: number;
  /** On-site damping Γ */
  gamma?: number;
  /** Preferred handedness for stability */
  preferredHandedness?: Handedness;
  /** Enable topological protection */
  topologicalProtection?: boolean;
  /** Enable CISS coherence enhancement */
  enableCISS?: boolean;
}

export interface ChiralState {
  /** Current handedness of dominant mode */
  handedness: Handedness;
  /** Chiral velocity c = η/Γ */
  chiralVelocity: number;
  /** Stability regime */
  stabilityRegime: 'stable' | 'transitional' | 'unstable';
  /** Topological winding number */
  windingNumber: number;
  /** Spin polarization (if CISS enabled) */
  spinPolarization: number;
  /** Non-reciprocal asymmetry measure */
  asymmetry: number;
  /** Protected mode amplitude */
  protectedAmplitude: number;
}

export interface ChiralCoupling {
  /** Forward coupling J_{i,j} */
  forward: number;
  /** Backward coupling J_{j,i} */
  backward: number;
  /** Asymmetry ratio */
  ratio: number;
}

// ============================================
// CHIRAL DYNAMICS ENGINE
// ============================================

export class ChiralEngine {
  private eta: number;
  private gamma: number;
  private preferredHandedness: Handedness;
  private topologicalProtection: boolean;
  private enableCISS: boolean;

  // Internal state
  private leftModeAmplitude: number = 0;
  private rightModeAmplitude: number = 0;
  private windingNumber: number = 0;
  private spinState: Float64Array;

  constructor(config: ChiralConfig = {}) {
    this.eta = config.eta ?? CHIRAL_ETA_DEFAULT;
    this.gamma = config.gamma ?? 1.0;
    this.preferredHandedness = config.preferredHandedness ?? 'right';
    this.topologicalProtection = config.topologicalProtection ?? true;
    this.enableCISS = config.enableCISS ?? true;
    this.spinState = new Float64Array(2); // [up, down] spin components
  }

  /**
   * Calculate non-reciprocal coupling between two indices
   * Implements J_{i,j} ≠ J_{j,i} for directional energy flow
   */
  calculateCoupling(i: number, j: number, baseStrength: number): ChiralCoupling {
    const direction = Math.sign(j - i);
    const distance = Math.abs(j - i);

    // Base coupling decays with distance
    const baseCoupling = baseStrength * Math.exp(-distance / PHI);

    // Non-reciprocal asymmetry based on chirality
    const asymmetry = this.eta * NON_RECIPROCAL_ASYMMETRY;

    // Forward and backward couplings differ
    const forward = baseCoupling * (1 + direction * asymmetry);
    const backward = baseCoupling * (1 - direction * asymmetry);

    return {
      forward,
      backward,
      ratio: forward / (backward + 1e-10),
    };
  }

  /**
   * Apply chiral dynamics to a state vector
   * Implements: φₜₜ - φₓₓ + sin(φ) = -ηφₓ - Γφₜ
   */
  applyChiralDynamics(
    state: Float64Array,
    velocity: Float64Array,
    dt: number
  ): { state: Float64Array; velocity: Float64Array } {
    const n = state.length;
    const newState = new Float64Array(n);
    const newVelocity = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      // Spatial derivative (central difference)
      const prev = state[(i - 1 + n) % n];
      const next = state[(i + 1) % n];
      const laplacian = prev - 2 * state[i] + next;

      // First spatial derivative for non-reciprocal term
      const gradient = (next - prev) / 2;

      // Sine-Gordon nonlinearity
      const nonlinear = Math.sin(state[i]);

      // Non-reciprocal dissipative dynamics
      // φₜₜ = φₓₓ - sin(φ) - ηφₓ - Γφₜ
      const acceleration =
        laplacian -
        nonlinear -
        this.eta * gradient -    // Non-reciprocal term
        this.gamma * velocity[i]; // Damping term

      // Update velocity and state
      newVelocity[i] = velocity[i] + acceleration * dt;
      newState[i] = state[i] + newVelocity[i] * dt;
    }

    // Apply topological protection if enabled
    if (this.topologicalProtection) {
      this.applyTopologicalProtection(newState, newVelocity);
    }

    // Update internal chiral mode tracking
    this.updateChiralModes(newState);

    return { state: newState, velocity: newVelocity };
  }

  /**
   * Apply topological protection to edge modes
   * Chiral edge states are robust to disorder
   */
  private applyTopologicalProtection(
    state: Float64Array,
    velocity: Float64Array
  ): void {
    const n = state.length;
    const edgeWidth = Math.ceil(n * 0.1); // 10% edge region

    // Protect edge modes from excessive perturbation
    for (let i = 0; i < edgeWidth; i++) {
      const protection = TOPOLOGICAL_PROTECTION_FACTOR *
        (1 - i / edgeWidth);

      // Left edge
      state[i] *= (1 - protection) + protection * Math.sign(state[i]);
      velocity[i] *= (1 - protection * 0.5);

      // Right edge
      const ri = n - 1 - i;
      state[ri] *= (1 - protection) + protection * Math.sign(state[ri]);
      velocity[ri] *= (1 - protection * 0.5);
    }
  }

  /**
   * Update tracking of left/right chiral modes
   */
  private updateChiralModes(state: Float64Array): void {
    const n = state.length;
    let leftSum = 0;
    let rightSum = 0;

    // Decompose into left and right propagating components
    for (let k = 1; k < n / 2; k++) {
      let realPart = 0;
      let imagPart = 0;

      for (let i = 0; i < n; i++) {
        const angle = (TAU * k * i) / n;
        realPart += state[i] * Math.cos(angle);
        imagPart += state[i] * Math.sin(angle);
      }

      const amplitude = Math.sqrt(realPart * realPart + imagPart * imagPart) / n;
      const phase = Math.atan2(imagPart, realPart);

      // Classify by phase velocity direction
      if (phase > 0) {
        rightSum += amplitude;
      } else {
        leftSum += amplitude;
      }
    }

    this.leftModeAmplitude = leftSum;
    this.rightModeAmplitude = rightSum;

    // Calculate winding number (topological invariant)
    this.windingNumber = this.calculateWindingNumber(state);
  }

  /**
   * Calculate topological winding number
   * Counts how many times the phase wraps around
   */
  private calculateWindingNumber(state: Float64Array): number {
    const n = state.length;
    let totalPhaseChange = 0;

    for (let i = 0; i < n; i++) {
      const current = state[i];
      const next = state[(i + 1) % n];

      // Phase difference
      let dPhase = Math.atan2(Math.sin(next - current), Math.cos(next - current));
      totalPhaseChange += dPhase;
    }

    return Math.round(totalPhaseChange / TAU);
  }

  /**
   * Apply CISS (Chiral-Induced Spin Selectivity) enhancement
   * Enhances coherence through spin-selective transport
   */
  applyCISSEnhancement(coherence: number): number {
    if (!this.enableCISS) return coherence;

    // Spin polarization based on handedness asymmetry
    const asymmetry = Math.abs(this.rightModeAmplitude - this.leftModeAmplitude) /
      (this.rightModeAmplitude + this.leftModeAmplitude + 1e-10);

    // Update spin state
    this.spinState[0] = 0.5 + asymmetry * CISS_COUPLING.maxPolarization / 2;
    this.spinState[1] = 0.5 - asymmetry * CISS_COUPLING.maxPolarization / 2;

    // Coherence boost from spin selectivity
    const boost = 1 + (CISS_COUPLING.coherenceBoost - 1) * asymmetry;

    // Suppress spin-flip decoherence
    const flipSuppression = 1 - CISS_COUPLING.flipSuppression * (1 - asymmetry);

    return coherence * boost * flipSuppression;
  }

  /**
   * Get current chiral velocity c = η/Γ
   * Determines propagation direction and stability regime
   */
  getChiralVelocity(): number {
    return this.eta / (this.gamma + 1e-10);
  }

  /**
   * Get current stability regime based on chiral velocity
   */
  getStabilityRegime(): 'stable' | 'transitional' | 'unstable' {
    const absVelocity = Math.abs(this.getChiralVelocity());

    if (absVelocity <= CHIRAL_STABILITY.stable.maxRatio) {
      return 'stable';
    } else if (absVelocity <= CHIRAL_STABILITY.transitional.maxRatio) {
      return 'transitional';
    }
    return 'unstable';
  }

  /**
   * Get current dominant handedness
   */
  getHandedness(): Handedness {
    const diff = this.rightModeAmplitude - this.leftModeAmplitude;
    const threshold = 0.1 * (this.rightModeAmplitude + this.leftModeAmplitude);

    if (Math.abs(diff) < threshold) return 'achiral';
    return diff > 0 ? 'right' : 'left';
  }

  /**
   * Get full chiral state
   */
  getState(): ChiralState {
    const velocity = this.getChiralVelocity();
    const totalAmplitude = this.leftModeAmplitude + this.rightModeAmplitude;

    return {
      handedness: this.getHandedness(),
      chiralVelocity: velocity,
      stabilityRegime: this.getStabilityRegime(),
      windingNumber: this.windingNumber,
      spinPolarization: this.spinState[0] - this.spinState[1],
      asymmetry: Math.abs(this.rightModeAmplitude - this.leftModeAmplitude) /
        (totalAmplitude + 1e-10),
      protectedAmplitude: this.preferredHandedness === 'right'
        ? this.rightModeAmplitude
        : this.leftModeAmplitude,
    };
  }

  /**
   * Adjust non-reciprocity to favor stability
   * Based on: waves with c = η/Γ are stable
   */
  optimizeForStability(targetVelocity: number = 0.618): void {
    // Adjust η to achieve target chiral velocity
    this.eta = targetVelocity * this.gamma;

    // Keep within golden ratio bounds for harmonic stability
    this.eta = Math.max(PHI_INVERSE * 0.5, Math.min(PHI, this.eta));
  }

  /**
   * Apply chiral gain to selectively amplify preferred handedness
   * Based on chiral gain-induced symmetry breaking research
   */
  applyChiralGain(
    state: Float64Array,
    gainStrength: number = 0.1
  ): Float64Array {
    const n = state.length;
    const amplified = new Float64Array(n);

    // Decompose, amplify preferred chirality, recompose
    for (let i = 0; i < n; i++) {
      // Local chirality measure (gradient direction)
      const prev = state[(i - 1 + n) % n];
      const next = state[(i + 1) % n];
      const localChirality = next - prev; // Positive = right-moving

      // Selective amplification based on preferred handedness
      const gain = this.preferredHandedness === 'right'
        ? gainStrength * Math.max(0, localChirality)
        : gainStrength * Math.max(0, -localChirality);

      amplified[i] = state[i] * (1 + gain);
    }

    return amplified;
  }

  /**
   * Calculate non-Hermitian effective Hamiltonian coupling
   * H_eff = H_S + Σ(J_{i,j} - iγ_{i,j}/2)O_i†O_j
   */
  getNonHermitianCoupling(i: number, j: number): { real: number; imag: number } {
    const baseCoupling = this.calculateCoupling(i, j, 1.0);

    // Real part: exchange coupling J
    const realPart = baseCoupling.forward;

    // Imaginary part: -iγ/2 for decay
    const imagPart = -NON_HERMITIAN.gammaCoupling / 2 *
      Math.cos(NON_HERMITIAN.phaseAngle * (j - i));

    return { real: realPart, imag: imagPart };
  }

  /**
   * Reset chiral engine state
   */
  reset(): void {
    this.leftModeAmplitude = 0;
    this.rightModeAmplitude = 0;
    this.windingNumber = 0;
    this.spinState.fill(0.5);
  }
}

// ============================================
// CONVENIENCE EXPORTS
// ============================================

export const createChiralEngine = (config?: ChiralConfig) =>
  new ChiralEngine(config);

/**
 * Quick check if a velocity ratio is in stable regime
 */
export const isChirallyStable = (eta: number, gamma: number): boolean => {
  const ratio = Math.abs(eta / (gamma + 1e-10));
  return ratio <= CHIRAL_VELOCITY_THRESHOLD;
};

/**
 * Calculate optimal eta for given gamma to achieve stability
 */
export const optimalEta = (gamma: number): number => {
  return PHI_INVERSE * gamma; // Golden ratio inverse for harmonic chirality
};
