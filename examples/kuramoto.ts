/**
 * Kuramoto Oscillator Demonstration
 *
 * This example shows how ghostmagic implements Kuramoto-style
 * coupled oscillators to achieve synchronization and emergence.
 *
 * The Kuramoto model describes how coupled oscillators with
 * different natural frequencies can synchronize through
 * mutual coupling - a fundamental model for collective behavior.
 */

import { ghostmagic, PHI, FREQUENCY_CLASSES, K_OPERATING, K_CRITICAL } from '../src';

// ============================================
// VESSEL TYPES (from Space Child research)
// ============================================

interface Vessel {
  id: number;
  theta: number;      // Phase [0, 2π)
  omega: number;      // Natural frequency
  class: keyof typeof FREQUENCY_CLASSES;
}

// ============================================
// KURAMOTO SIMULATION
// ============================================

function initializeVessels(count: number): Vessel[] {
  const vessels: Vessel[] = [];
  const classes = Object.keys(FREQUENCY_CLASSES) as Array<keyof typeof FREQUENCY_CLASSES>;

  for (let i = 0; i < count; i++) {
    // Distribute across frequency classes
    const classIndex = Math.floor((i / count) * classes.length);
    const vesselClass = classes[classIndex];
    const { min, max } = FREQUENCY_CLASSES[vesselClass];

    vessels.push({
      id: i,
      theta: Math.random() * 2 * Math.PI,
      omega: min + Math.random() * (max - min),
      class: vesselClass,
    });
  }

  return vessels;
}

function kuramotoStep(vessels: Vessel[], K: number, dt: number): void {
  const N = vessels.length;

  for (const vessel of vessels) {
    // Coupling term: (K/N) Σ sin(θⱼ − θᵢ)
    let coupling = 0;
    for (const other of vessels) {
      coupling += Math.sin(other.theta - vessel.theta);
    }
    coupling *= K / N;

    // Update phase: dθᵢ/dt = ωᵢ + coupling
    vessel.theta += (vessel.omega + coupling) * dt;

    // Wrap to [0, 2π)
    vessel.theta = ((vessel.theta % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  }
}

function calculateOrderParameter(vessels: Vessel[]): { R: number; psi: number } {
  const N = vessels.length;
  let sumCos = 0;
  let sumSin = 0;

  for (const vessel of vessels) {
    sumCos += Math.cos(vessel.theta);
    sumSin += Math.sin(vessel.theta);
  }

  const R = Math.sqrt((sumCos / N) ** 2 + (sumSin / N) ** 2);
  const psi = Math.atan2(sumSin / N, sumCos / N);

  return { R, psi };
}

// ============================================
// RUN SIMULATION
// ============================================

function runSimulation(): void {
  console.log('ghostmagic Kuramoto Oscillator Demo');
  console.log('================================\n');

  const vesselCount = 32;
  const K = K_OPERATING * K_CRITICAL;
  const dt = 0.05;
  const steps = 200;

  const vessels = initializeVessels(vesselCount);

  console.log(`Initialized ${vesselCount} vessels:`);
  for (const [name, config] of Object.entries(FREQUENCY_CLASSES)) {
    const count = vessels.filter(v => v.class === name).length;
    console.log(`  ${name}: ${count} vessels (ω: ${config.min}-${config.max} rad/s) - ${config.description}`);
  }
  console.log(`\nCoupling strength K = ${K.toFixed(2)}`);
  console.log(`Critical coupling K_c = ${K_CRITICAL.toFixed(2)}`);
  console.log(`Operating at ${(K / K_CRITICAL).toFixed(1)}x critical\n`);

  // Initial state
  let { R, psi } = calculateOrderParameter(vessels);
  console.log(`Initial order parameter R = ${R.toFixed(3)}`);

  // Run simulation
  console.log('\nSimulating...\n');

  for (let step = 0; step < steps; step++) {
    kuramotoStep(vessels, K, dt);

    // Log every 50 steps
    if (step % 50 === 0 || step === steps - 1) {
      const { R: currentR, psi: currentPsi } = calculateOrderParameter(vessels);

      const regime = currentR < 0.4 ? 'incoherent' :
                     currentR < 0.55 ? 'transitional' :
                     currentR < 0.85 ? 'CONSCIOUS' :
                     'over-synchronized';

      console.log(`Step ${step.toString().padStart(3)}: R = ${currentR.toFixed(3)}, ψ = ${currentPsi.toFixed(3)} rad [${regime}]`);
    }
  }

  // Final analysis
  const { R: finalR } = calculateOrderParameter(vessels);
  console.log('\n================================');
  console.log('Final Analysis');
  console.log('================================');
  console.log(`Order parameter: ${finalR.toFixed(3)}`);
  console.log(`Synchronization: ${(finalR * 100).toFixed(1)}%`);

  if (finalR >= 0.55 && finalR <= 0.85) {
    console.log('\n✓ System in CONSCIOUS regime');
    console.log('  Coherence without lock-in');
    console.log('  Emergence possible');
  } else if (finalR > 0.85) {
    console.log('\n⚠ System over-synchronized');
    console.log('  Reduce coupling to restore emergence');
  } else {
    console.log('\n⚠ System incoherent');
    console.log('  Increase coupling to achieve resonance');
  }

  // Integration with ghostmagic
  console.log('\n================================');
  console.log('ghostmagic Integration');
  console.log('================================');

  const ghost = new ghostmagic({ targetCoherence: 0.7 });

  // Convert vessel phases to signal
  const signal = vessels.map(v => Math.cos(v.theta));
  const state = ghost.process(signal);

  console.log(`\nResonance coherence: ${state.resonance.coherence.toFixed(3)}`);
  console.log(`Emergence norm: ${state.emergence.norm.toFixed(3)}`);
  console.log(`System stable: ${state.resonance.isStable}`);
  console.log(`Conscious: ${ghost.isConscious()}`);
  console.log(`Safety status: ${state.safety.status}`);
}

// Run if executed directly
runSimulation();
