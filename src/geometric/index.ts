/**
 * ghostOS Geometric Control Layer
 *
 * Differential geometry primitives for the oscillator phase space.
 * Provides Berry curvature, manifold operations, curvature tensors,
 * and chiral anomaly detection for the resonant systems architecture.
 *
 * Module overview:
 * - manifold: OscillatorManifold class with BFGS metric adaptation on T^N
 * - berry: Berry phase, connection, and curvature computation
 * - curvature: Full differential geometry pipeline (Christoffel -> Riemann -> Ricci)
 * - anomaly: Chiral anomaly detection and spectral asymmetry monitoring
 */

// Manifold operations
export {
  OscillatorManifold,
  type ManifoldPoint,
  type ManifoldMetrics,
} from './manifold';

// Berry phase and connection
export {
  computeBerryConnection,
  computeBerryCurvature,
  computeBerryPhaseLoop,
  createBerryTracker,
  type BerryState,
} from './berry';

// Curvature tensors
export {
  computeChristoffelSymbols,
  computeRiemannTensor,
  computeRicciTensor,
  computeRicciScalar,
  invertMetric,
  computeMetricDerivatives,
  computeScalarCurvature,
} from './curvature';

// Chiral anomaly detection
export {
  detectChiralAnomaly,
  computeSpectralAsymmetry,
  computeTopologicalCharge,
  createAnomalyMonitor,
  type ChiralAnomalyResult,
} from './anomaly';
