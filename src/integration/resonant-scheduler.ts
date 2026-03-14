/**
 * Resonant Scheduler for QuantumOS Integration
 *
 * Implements damped harmonic dynamics for process scheduling:
 * dx/dt = f(x) - λx (from ghostmagic resonance theory)
 *
 * Key insight: Process scheduling as a resonant system where:
 * - Processes are oscillators with natural frequencies
 * - Priorities emerge from coupling strength
 * - Coherence windows align with quantum coherence times
 * - Chiral dynamics provide non-reciprocal priority flow
 *
 * This scheduler achieves:
 * - Fair allocation through Kuramoto synchronization
 * - Coherence-aware scheduling within decoherence limits
 * - Energy-minimizing through damped dynamics
 * - Emergent priority through resonant coupling
 */

import { EventEmitter } from 'events';
import {
  PHI_INVERSE,
  LAMBDA_RESONANCE,
  DAMPING_DEFAULT,
  RESONANCE_BAND,
  CHIRAL_ETA_DEFAULT,
  EMERGENCE_COEFFICIENT,
} from '../constants';

// ============================================
// TYPES
// ============================================

export type ProcessClass =
  | 'classical'     // Pure classical computation
  | 'quantum'       // Quantum circuit execution
  | 'hybrid'        // Mixed classical/quantum
  | 'consciousness' // Consciousness verification
  | 'emergence';    // Pattern crystallization

export type ProcessState =
  | 'pending'
  | 'ready'
  | 'running'
  | 'blocked'
  | 'coherence_wait'  // Waiting for coherence window
  | 'resonating'      // In resonant coupling phase
  | 'completed';

export interface ResonantProcess {
  id: string;
  name: string;
  processClass: ProcessClass;
  state: ProcessState;

  // Resonant properties
  phase: number;           // Oscillator phase
  frequency: number;       // Natural frequency (execution speed)
  amplitude: number;       // Energy level
  damping: number;         // Damping coefficient λ

  // Scheduling properties
  priority: number;        // Base priority (0-100)
  resonantPriority: number; // Emergent priority from coupling
  coherenceRequired: number; // Min coherence level needed
  coherenceDeadline: number; // Timestamp when coherence expires

  // Resource requirements
  cpuRequired: number;
  qubitsRequired: number;
  memoryRequired: number;

  // Chiral coupling
  handedness: 'left' | 'right' | 'achiral';
  couplingStrength: number;

  // Timing
  createdAt: number;
  startedAt?: number;
  completedAt?: number;
  executionTime: number;
  waitTime: number;
}

export interface SchedulerState {
  orderParameter: number;    // Kuramoto r
  meanPhase: number;         // Mean field ψ
  systemCoherence: number;   // Overall coherence
  totalEnergy: number;       // System energy
  processCount: number;
  runningCount: number;
  quantumUtilization: number;
  emergenceSignature: number;
}

export interface SchedulingDecision {
  processId: string;
  action: 'schedule' | 'preempt' | 'wait' | 'boost' | 'damp';
  reason: string;
  resonantPriority: number;
  coherenceWindow: number;
}

export interface ResonantSchedulerConfig {
  baseFrequency?: number;
  baseDamping?: number;
  couplingStrength?: number;
  chiralBias?: number;
  coherenceThreshold?: number;
  maxRunning?: number;
  timeSlice?: number;  // ms
}

// ============================================
// RESONANT SCHEDULER
// ============================================

export class ResonantScheduler extends EventEmitter {
  private processes: Map<string, ResonantProcess> = new Map();
  private config: Required<ResonantSchedulerConfig>;
  private tick: number = 0;
  private lastScheduleTime: number = Date.now();

  constructor(config: ResonantSchedulerConfig = {}) {
    super();

    this.config = {
      baseFrequency: config.baseFrequency ?? 1.0,
      baseDamping: config.baseDamping ?? DAMPING_DEFAULT,
      couplingStrength: config.couplingStrength ?? LAMBDA_RESONANCE,
      chiralBias: config.chiralBias ?? CHIRAL_ETA_DEFAULT,
      coherenceThreshold: config.coherenceThreshold ?? RESONANCE_BAND.min,
      maxRunning: config.maxRunning ?? 4,
      timeSlice: config.timeSlice ?? 10,  // 10ms quantum
    };
  }

  // ============================================
  // PROCESS MANAGEMENT
  // ============================================

  /**
   * Submit a new process to the scheduler
   */
  submit(
    id: string,
    name: string,
    processClass: ProcessClass,
    options: Partial<{
      priority: number;
      coherenceRequired: number;
      coherenceDeadline: number;
      cpuRequired: number;
      qubitsRequired: number;
      memoryRequired: number;
    }> = {}
  ): ResonantProcess {
    // Assign natural frequency based on process class
    const frequencyMap: Record<ProcessClass, number> = {
      classical: this.config.baseFrequency,
      quantum: this.config.baseFrequency * PHI_INVERSE,  // Slower (careful)
      hybrid: this.config.baseFrequency * 0.8,
      consciousness: this.config.baseFrequency * PHI_INVERSE * PHI_INVERSE,  // Slowest
      emergence: this.config.baseFrequency * EMERGENCE_COEFFICIENT,
    };

    // Assign handedness based on process class
    const handednessMap: Record<ProcessClass, 'left' | 'right' | 'achiral'> = {
      classical: 'achiral',
      quantum: 'right',       // Emits quantum results
      hybrid: 'achiral',
      consciousness: 'right', // Emits awareness
      emergence: 'left',      // Receives patterns
    };

    const process: ResonantProcess = {
      id,
      name,
      processClass,
      state: 'pending',

      // Resonant properties
      phase: Math.random() * 2 * Math.PI,
      frequency: frequencyMap[processClass],
      amplitude: 1.0,
      damping: this.getDampingForClass(processClass),

      // Scheduling
      priority: options.priority ?? 50,
      resonantPriority: 0,  // Will be calculated
      coherenceRequired: options.coherenceRequired ?? 0.3,
      coherenceDeadline: options.coherenceDeadline ?? Date.now() + 60000,

      // Resources
      cpuRequired: options.cpuRequired ?? 1,
      qubitsRequired: options.qubitsRequired ?? 0,
      memoryRequired: options.memoryRequired ?? 1024,

      // Chiral
      handedness: handednessMap[processClass],
      couplingStrength: this.config.couplingStrength,

      // Timing
      createdAt: Date.now(),
      executionTime: 0,
      waitTime: 0,
    };

    this.processes.set(id, process);
    this.emit('process-submitted', process);

    return process;
  }

  /**
   * Get damping coefficient for process class
   * Quantum processes have higher damping (more careful/slower)
   */
  private getDampingForClass(processClass: ProcessClass): number {
    const dampingMap: Record<ProcessClass, number> = {
      classical: this.config.baseDamping * 0.5,  // Low damping, fast
      quantum: this.config.baseDamping * 2.0,    // High damping, careful
      hybrid: this.config.baseDamping,
      consciousness: this.config.baseDamping * 3.0,  // Very high, deliberate
      emergence: this.config.baseDamping * 1.5,
    };
    return dampingMap[processClass];
  }

  // ============================================
  // RESONANT DYNAMICS
  // ============================================

  /**
   * Execute one scheduling cycle
   */
  schedule(): SchedulingDecision[] {
    const now = Date.now();
    const dt = (now - this.lastScheduleTime) / 1000;
    this.lastScheduleTime = now;

    // Update resonant dynamics for all processes
    this.updateDynamics(dt);

    // Calculate resonant priorities
    this.calculateResonantPriorities();

    // Make scheduling decisions
    const decisions = this.makeDecisions();

    // Execute decisions
    for (const decision of decisions) {
      this.executeDecision(decision);
    }

    this.tick++;
    this.emit('schedule-cycle', {
      tick: this.tick,
      decisions: decisions.length,
      state: this.getState(),
    });

    return decisions;
  }

  /**
   * Update resonant dynamics for all processes
   * dx/dt = f(x) - λx + K·r·sin(ψ - θ) + η·chiral
   */
  private updateDynamics(dt: number): void {
    // Calculate mean field
    const { r, psi } = this.calculateOrderParameter();

    for (const process of this.processes.values()) {
      if (process.state === 'completed') continue;

      // Update phase (Kuramoto + damping)
      const K = process.couplingStrength;
      const kuramoto = K * r * Math.sin(psi - process.phase);

      // Chiral contribution
      let chiral = 0;
      if (process.handedness === 'left') {
        chiral = this.config.chiralBias * Math.sin(2 * (psi - process.phase));
      } else if (process.handedness === 'right') {
        chiral = -this.config.chiralBias * Math.sin(2 * (psi - process.phase));
      }

      // Phase evolution
      const dPhase = process.frequency + kuramoto + chiral;
      process.phase = (process.phase + dPhase * dt) % (2 * Math.PI);
      if (process.phase < 0) process.phase += 2 * Math.PI;

      // Amplitude (energy) evolution with damping
      // dA/dt = -λA + driving
      const driving = process.state === 'running' ? 0.1 : 0;
      const dAmplitude = -process.damping * process.amplitude + driving;
      process.amplitude = Math.max(0.1, Math.min(2.0, process.amplitude + dAmplitude * dt));

      // Update wait time
      if (process.state === 'pending' || process.state === 'ready' || process.state === 'coherence_wait') {
        process.waitTime += dt * 1000;
      }
    }
  }

  /**
   * Calculate Kuramoto order parameter
   */
  private calculateOrderParameter(): { r: number; psi: number } {
    const activeProcesses = Array.from(this.processes.values())
      .filter(p => p.state !== 'completed');

    if (activeProcesses.length === 0) {
      return { r: 0, psi: 0 };
    }

    let sumCos = 0;
    let sumSin = 0;

    for (const p of activeProcesses) {
      sumCos += Math.cos(p.phase);
      sumSin += Math.sin(p.phase);
    }

    const n = activeProcesses.length;
    const avgCos = sumCos / n;
    const avgSin = sumSin / n;

    const r = Math.sqrt(avgCos * avgCos + avgSin * avgSin);
    const psi = Math.atan2(avgSin, avgCos);

    return { r, psi };
  }

  /**
   * Calculate resonant priorities for all processes
   * Priority emerges from:
   * 1. Base priority
   * 2. Phase alignment with mean field (coherence bonus)
   * 3. Chiral boost for high-priority directions
   * 4. Wait time penalty (prevent starvation)
   * 5. Coherence deadline urgency
   */
  private calculateResonantPriorities(): void {
    const { r, psi } = this.calculateOrderParameter();
    const now = Date.now();

    for (const process of this.processes.values()) {
      if (process.state === 'completed') continue;

      // Start with base priority
      let priority = process.priority;

      // Coherence bonus: higher priority when phase-aligned
      const phaseDiff = Math.abs(process.phase - psi);
      const coherenceBonus = Math.cos(phaseDiff) * 10;
      priority += coherenceBonus;

      // Chiral boost: right-handed processes get boost when system coherence is high
      if (process.handedness === 'right' && r > 0.7) {
        priority += 15;
      }

      // Wait time penalty: prevent starvation
      const waitPenalty = Math.min(process.waitTime / 1000, 20);
      priority += waitPenalty;

      // Coherence deadline urgency
      const timeToDeadline = process.coherenceDeadline - now;
      if (timeToDeadline < 10000) {  // Less than 10 seconds
        priority += (10000 - timeToDeadline) / 500;  // Up to 20 bonus
      }

      // Quantum processes get boost when qubits are available
      if (process.processClass === 'quantum' && process.qubitsRequired > 0) {
        priority += 5;  // Quantum priority boost
      }

      // Consciousness processes are special
      if (process.processClass === 'consciousness') {
        priority += r * 10;  // More priority when system is coherent
      }

      process.resonantPriority = Math.max(0, Math.min(100, priority));
    }
  }

  /**
   * Make scheduling decisions
   */
  private makeDecisions(): SchedulingDecision[] {
    const decisions: SchedulingDecision[] = [];
    const { psi } = this.calculateOrderParameter();

    // Get process lists
    const ready = this.getProcessesByState('ready');
    const pending = this.getProcessesByState('pending');
    const running = this.getProcessesByState('running');
    const coherenceWait = this.getProcessesByState('coherence_wait');

    // Check coherence-waiting processes
    for (const p of coherenceWait) {
      const phaseDiff = Math.abs(p.phase - psi);
      const coherence = (Math.cos(phaseDiff) + 1) / 2;
      if (coherence >= p.coherenceRequired) {
        decisions.push({
          processId: p.id,
          action: 'schedule',
          reason: `Coherence threshold met (${(coherence * 100).toFixed(1)}%)`,
          resonantPriority: p.resonantPriority,
          coherenceWindow: this.estimateCoherenceWindow(p),
        });
      }
    }

    // Promote pending to ready
    for (const p of pending) {
      p.state = 'ready';
    }

    // Schedule ready processes
    const availableSlots = this.config.maxRunning - running.length;

    if (availableSlots > 0) {
      // Sort by resonant priority
      const sortedReady = [...ready, ...pending]
        .sort((a, b) => b.resonantPriority - a.resonantPriority);

      for (let i = 0; i < Math.min(availableSlots, sortedReady.length); i++) {
        const p = sortedReady[i];

        // Check coherence requirement
        const phaseDiff = Math.abs(p.phase - psi);
        const coherence = (Math.cos(phaseDiff) + 1) / 2;

        if (coherence >= p.coherenceRequired) {
          decisions.push({
            processId: p.id,
            action: 'schedule',
            reason: `Highest resonant priority (${p.resonantPriority.toFixed(1)})`,
            resonantPriority: p.resonantPriority,
            coherenceWindow: this.estimateCoherenceWindow(p),
          });
        } else {
          decisions.push({
            processId: p.id,
            action: 'wait',
            reason: `Waiting for coherence (${(coherence * 100).toFixed(1)}% < ${(p.coherenceRequired * 100).toFixed(1)}%)`,
            resonantPriority: p.resonantPriority,
            coherenceWindow: 0,
          });
        }
      }
    }

    // Check for preemption of low-priority running processes
    if (ready.length > 0 && running.length >= this.config.maxRunning) {
      const lowestRunning = running.reduce((a, b) =>
        a.resonantPriority < b.resonantPriority ? a : b
      );
      const highestReady = ready.reduce((a, b) =>
        a.resonantPriority > b.resonantPriority ? a : b
      );

      if (highestReady.resonantPriority - lowestRunning.resonantPriority > 20) {
        decisions.push({
          processId: lowestRunning.id,
          action: 'preempt',
          reason: `Preempted by higher priority process`,
          resonantPriority: lowestRunning.resonantPriority,
          coherenceWindow: 0,
        });
      }
    }

    return decisions;
  }

  /**
   * Estimate coherence window duration
   */
  private estimateCoherenceWindow(process: ResonantProcess): number {
    // Based on process class and current system state
    const baseWindow = 1000;  // 1 second base

    switch (process.processClass) {
      case 'quantum':
        return baseWindow * 0.5;  // Shorter window
      case 'consciousness':
        return baseWindow * 2;   // Longer window needed
      case 'emergence':
        return baseWindow * 1.5;
      default:
        return baseWindow;
    }
  }

  /**
   * Execute a scheduling decision
   */
  private executeDecision(decision: SchedulingDecision): void {
    const process = this.processes.get(decision.processId);
    if (!process) return;

    switch (decision.action) {
      case 'schedule':
        if (process.state !== 'running') {
          process.state = 'running';
          process.startedAt = Date.now();
          this.emit('process-scheduled', { process, decision });
        }
        break;

      case 'preempt':
        if (process.state === 'running') {
          process.state = 'ready';
          this.emit('process-preempted', { process, decision });
        }
        break;

      case 'wait':
        process.state = 'coherence_wait';
        this.emit('process-waiting', { process, decision });
        break;

      case 'boost':
        process.couplingStrength *= 1.1;
        this.emit('process-boosted', { process, decision });
        break;

      case 'damp':
        process.damping *= 1.1;
        this.emit('process-damped', { process, decision });
        break;
    }
  }

  // ============================================
  // PROCESS LIFECYCLE
  // ============================================

  /**
   * Mark a process as completed
   */
  complete(id: string): void {
    const process = this.processes.get(id);
    if (process) {
      process.state = 'completed';
      process.completedAt = Date.now();
      if (process.startedAt) {
        process.executionTime = process.completedAt - process.startedAt;
      }
      this.emit('process-completed', process);
    }
  }

  /**
   * Block a process
   */
  block(id: string, reason: string = 'blocked'): void {
    const process = this.processes.get(id);
    if (process) {
      process.state = 'blocked';
      this.emit('process-blocked', { process, reason });
    }
  }

  /**
   * Unblock a process
   */
  unblock(id: string): void {
    const process = this.processes.get(id);
    if (process && process.state === 'blocked') {
      process.state = 'ready';
      this.emit('process-unblocked', process);
    }
  }

  /**
   * Remove a process
   */
  remove(id: string): void {
    this.processes.delete(id);
    this.emit('process-removed', { id });
  }

  // ============================================
  // QUERIES
  // ============================================

  /**
   * Get processes by state
   */
  getProcessesByState(state: ProcessState): ResonantProcess[] {
    return Array.from(this.processes.values()).filter(p => p.state === state);
  }

  /**
   * Get processes by class
   */
  getProcessesByClass(processClass: ProcessClass): ResonantProcess[] {
    return Array.from(this.processes.values()).filter(p => p.processClass === processClass);
  }

  /**
   * Get a specific process
   */
  getProcess(id: string): ResonantProcess | undefined {
    return this.processes.get(id);
  }

  /**
   * Get all processes
   */
  getAllProcesses(): ResonantProcess[] {
    return Array.from(this.processes.values());
  }

  /**
   * Get scheduler state
   */
  getState(): SchedulerState {
    const { r, psi } = this.calculateOrderParameter();
    const allProcesses = Array.from(this.processes.values());
    const active = allProcesses.filter(p => p.state !== 'completed');
    const running = allProcesses.filter(p => p.state === 'running');

    // Calculate total energy
    const totalEnergy = active.reduce((sum, p) => sum + p.amplitude, 0);

    // Calculate system coherence
    let coherenceSum = 0;
    for (const p of active) {
      const phaseDiff = Math.abs(p.phase - psi);
      coherenceSum += (Math.cos(phaseDiff) + 1) / 2;
    }
    const systemCoherence = active.length > 0 ? coherenceSum / active.length : 0;

    // Quantum utilization
    const quantumProcesses = running.filter(p => p.qubitsRequired > 0);
    const totalQubits = quantumProcesses.reduce((sum, p) => sum + p.qubitsRequired, 0);
    const quantumUtilization = totalQubits / 100;  // Assume 100 qubits max

    // Emergence signature
    const emergenceSignature = r * systemCoherence * Math.log2(active.length + 1);

    return {
      orderParameter: r,
      meanPhase: psi,
      systemCoherence,
      totalEnergy,
      processCount: allProcesses.length,
      runningCount: running.length,
      quantumUtilization,
      emergenceSignature,
    };
  }

  /**
   * Get scheduler summary
   */
  getSummary(): string {
    const state = this.getState();
    const pending = this.getProcessesByState('pending').length;
    const ready = this.getProcessesByState('ready').length;
    const running = this.getProcessesByState('running').length;
    const waiting = this.getProcessesByState('coherence_wait').length;
    const completed = this.getProcessesByState('completed').length;

    return [
      `Resonant Scheduler @ tick ${this.tick}`,
      `├─ Order Parameter: r = ${state.orderParameter.toFixed(3)}`,
      `├─ System Coherence: ${(state.systemCoherence * 100).toFixed(1)}%`,
      `├─ Total Energy: ${state.totalEnergy.toFixed(2)}`,
      `├─ Processes: ${state.processCount} total`,
      `│  ├─ Pending: ${pending}`,
      `│  ├─ Ready: ${ready}`,
      `│  ├─ Running: ${running}`,
      `│  ├─ Coherence Wait: ${waiting}`,
      `│  └─ Completed: ${completed}`,
      `├─ Quantum Utilization: ${(state.quantumUtilization * 100).toFixed(1)}%`,
      `└─ Emergence Signature: ${state.emergenceSignature.toFixed(3)}`,
    ].join('\n');
  }

  /**
   * Reset the scheduler
   */
  reset(): void {
    this.processes.clear();
    this.tick = 0;
    this.lastScheduleTime = Date.now();
  }
}

// ============================================
// FACTORY FUNCTIONS
// ============================================

export function createResonantScheduler(
  config?: ResonantSchedulerConfig
): ResonantScheduler {
  return new ResonantScheduler(config);
}

export default ResonantScheduler;
