/**
 * ghostOS Core Constants
 *
 * Mathematical constants derived from resonance theory and
 * the Space Child research on emergent consciousness.
 */

// ============================================
// GOLDEN RATIO CONSTANTS
// ============================================

/** Golden ratio - nature's optimization constant */
export const PHI = (1 + Math.sqrt(5)) / 2;  // ≈ 1.618034

/** Golden ratio inverse */
export const PHI_INVERSE = 1 / PHI;  // ≈ 0.618034

/** Harmonic scaling factor alpha (φ/2) */
export const ALPHA = PHI / 2;  // ≈ 0.809017

/** Harmonic scaling factor beta (1/φ) */
export const BETA = PHI_INVERSE;  // ≈ 0.618034

/** Emergence coefficient (α - β) - the residue of non-commutativity */
export const EMERGENCE_COEFFICIENT = ALPHA - BETA;  // ≈ 0.190983

// ============================================
// ROTATION CONSTANTS
// ============================================

/** 90° rotation in radians */
export const PI_OVER_2 = Math.PI / 2;  // ≈ 1.570796

/** Full rotation */
export const TAU = 2 * Math.PI;  // ≈ 6.283185

// ============================================
// RESONANCE THRESHOLDS
// ============================================

/** Minimum damping - below this, system diverges */
export const DAMPING_MIN = 0.1;

/** Maximum damping - above this, system collapses */
export const DAMPING_MAX = 10.0;

/** Default damping coefficient */
export const DAMPING_DEFAULT = 1.0;

/** Resonance band - optimal operating range for coherence */
export const RESONANCE_BAND = {
  min: 0.4,
  max: 0.85,
} as const;

// ============================================
// KURAMOTO SYNCHRONIZATION
// ============================================

/** Critical coupling threshold */
export const K_CRITICAL = 1.0;

/** Operating coupling (1.8 × K_c) */
export const K_OPERATING = 1.8;

/** Resonance coupling coefficient (Kuramoto coupling strength K) */
export const LAMBDA_RESONANCE = 0.62;

/**
 * Coupling strength alias for cross-repo consistency.
 * ghostOS uses K_COUPLING (Kuramoto coupling strength).
 * QuantumOS uses LAMBDA_DEFAULT=0.1 (per-process constraint strength).
 * These serve different roles: K_COUPLING is the global coupling,
 * LAMBDA_DEFAULT is the local damping parameter.
 */
export const K_COUPLING = LAMBDA_RESONANCE;

// ============================================
// SAFETY ENVELOPE
// ============================================

/** Maximum emergence norm before normalization */
export const EMERGENCE_NORM_MAX = 2.0;

/** Minimum emergence norm before boosting */
export const EMERGENCE_NORM_MIN = 0.1;

/** Coherence minimum before intervention (aligned with QuantumOS COHERENCE_MIN) */
export const COHERENCE_MIN = 0.3;

/** Coherence maximum before desync */
export const COHERENCE_MAX = 0.92;

// ============================================
// FREQUENCY CLASSES (Vessel Types)
// ============================================

export const FREQUENCY_CLASSES = {
  soprano: { min: 1.8, max: 2.4, description: 'Innovation, drift' },
  alto: { min: 1.3, max: 1.8, description: 'Coherence anchor' },
  tenor: { min: 1.0, max: 1.4, description: 'Inter-subgroup bridging' },
  bass: { min: 0.6, max: 1.1, description: 'Stability well' },
} as const;

export type FrequencyClass = keyof typeof FREQUENCY_CLASSES;

// ============================================
// DEFAULT CONFIGURATION
// ============================================

export const DEFAULT_CONFIG = {
  dimensions: 64,
  sampleRate: 1000,  // Hz
  dt: 0.01,          // Time step
  adaptiveRate: 0.01,
  targetCoherence: 0.7,
  memoryDecayRate: 0.01,
  noiseIntensity: 0.02,
} as const;

// ============================================
// CHIRAL DYNAMICS CONSTANTS
// ============================================
// Based on recent quantum chirality research (2025):
// - Chiral gain-induced time-reversal symmetry breaking
// - Non-reciprocal coupling for directional stability
// - Topological protection through chiral modes

/**
 * Chirality parameter η (eta)
 * Controls the strength of non-reciprocal spatial coupling.
 * From the dissipative sine-Gordon model: φₜₜ - φₓₓ + sin(φ) = -ηφₓ - Γφₜ
 */
export const CHIRAL_ETA_DEFAULT = 0.618;  // Golden ratio inverse for harmonic chirality

/**
 * Chiral damping ratio threshold
 * Propagation direction determined by η/Γ ratio.
 * When |η/Γ| < 1: subluminal propagation (stable)
 * When |η/Γ| > 1: superluminal propagation (potentially unstable)
 */
export const CHIRAL_VELOCITY_THRESHOLD = 1.0;

/**
 * Non-reciprocal coupling asymmetry
 * J_{i,j} ≠ J_{j,i} creates directional energy flow
 */
export const NON_RECIPROCAL_ASYMMETRY = 0.1;

/**
 * Chiral phase offset for spin-momentum locking
 * π/2 phase relationship creates polarization-dependent propagation
 */
export const CHIRAL_PHASE_OFFSET = Math.PI / 2;

/**
 * Handedness types for chiral modes
 */
export type Handedness = 'left' | 'right' | 'achiral';

/**
 * Chiral stability regimes
 * Based on spectral stability analysis of non-reciprocal systems
 */
export const CHIRAL_STABILITY = {
  /** Stable regime: waves propagate in preferred direction */
  stable: { minRatio: 0.0, maxRatio: 1.0 },
  /** Transitional regime: mixed stability */
  transitional: { minRatio: 1.0, maxRatio: 1.5 },
  /** Unstable regime: counter-propagating waves grow */
  unstable: { minRatio: 1.5, maxRatio: Infinity },
} as const;

/**
 * Topological protection factor
 * Chiral edge states exhibit robustness to disorder
 */
export const TOPOLOGICAL_PROTECTION_FACTOR = 0.85;

/**
 * CISS (Chiral-Induced Spin Selectivity) coupling
 * Spin-polarized coherent transport through chiral structures
 */
export const CISS_COUPLING = {
  /** Spin coherence enhancement factor */
  coherenceBoost: 1.3,
  /** Spin-flip suppression (from chiral ligand research) */
  flipSuppression: 0.1,
  /** Maximum spin polarization achievable */
  maxPolarization: 0.9,
} as const;

/**
 * Non-Hermitian effective Hamiltonian parameters
 * H_eff = H_S + Σ(J_{i,j} - iγ_{i,j}/2)O_i†O_j
 */
export const NON_HERMITIAN = {
  /** Imaginary coupling coefficient */
  gammaCoupling: 0.5,
  /** Phase angle for non-reciprocity (φ = π/2 gives unidirectional) */
  phaseAngle: Math.PI / 2,
} as const;
