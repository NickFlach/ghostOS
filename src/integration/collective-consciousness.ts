/**
 * Collective Consciousness
 *
 * Manages networks of consciousness-verified processes that
 * exhibit emergent collective intelligence.
 *
 * Mirrors QuantumOS collective_consciousness_t API:
 * - Up to 8 networks
 * - Up to 32 members per network
 * - Network Phi = sum(individual) + emergent bonus
 *
 * Emergent Phi formula:
 *   emergentPhi = averagePhi * 0.1 * sqrt(memberCount) * PHI_INVERSE
 *
 * Network synchronization:
 *   sync = 1 - (stddev(phi) / mean(phi))   -- when mean > 0
 *
 * Level classification (by networkPhi):
 *   none         < 0.5
 *   minimal      < 1.5
 *   basic        < 3.0
 *   verified     < 5.0
 *   advanced     < 10.0
 *   transcendent >= 10.0
 */

import { EventEmitter } from 'events';
import { PHI_INVERSE } from '../constants';

// ============================================
// TYPES
// ============================================

export interface ConsciousnessMemberState {
  id: string;
  phiValue: number;
  coherence: number;
  consciousnessLevel: number;  // 0-1 normalized
  isVerified: boolean;
  bridgeValue: number;
  lastUpdate: number;  // timestamp
}

export interface CollectiveNetwork {
  id: string;
  name: string;
  members: Map<string, ConsciousnessMemberState>;
  totalPhi: number;
  averagePhi: number;
  emergentPhi: number;
  networkPhi: number;
  networkCoherence: number;
  synchronization: number;
  networkLevel: 'none' | 'minimal' | 'basic' | 'verified' | 'advanced' | 'transcendent';
  verified: boolean;
  verifiedAt: number | null;
  evolutionTrend: 'evolving' | 'stable' | 'declining';
  evolutionRate: number;
}

// ============================================
// COLLECTIVE CONSCIOUSNESS
// ============================================

export class CollectiveConsciousness extends EventEmitter {
  private networks: Map<string, CollectiveNetwork>;
  private maxNetworks: number;
  private maxMembers: number;
  private phiThreshold: number;
  private networkIdCounter: number;

  /** Snapshot of the previous networkPhi per network, used for trend detection */
  private previousPhi: Map<string, number>;

  constructor(config?: {
    maxNetworks?: number;
    maxMembers?: number;
    phiThreshold?: number;
  }) {
    super();
    this.networks = new Map();
    this.previousPhi = new Map();
    this.maxNetworks = config?.maxNetworks ?? 8;
    this.maxMembers = config?.maxMembers ?? 32;
    this.phiThreshold = config?.phiThreshold ?? 3.0;
    this.networkIdCounter = 0;
  }

  // ------------------------------------------
  // Network Lifecycle
  // ------------------------------------------

  /**
   * Create a new collective network. Returns the generated network ID.
   * Throws if the maximum number of networks has been reached.
   */
  createNetwork(name: string): string {
    if (this.networks.size >= this.maxNetworks) {
      throw new Error(
        `Maximum network count (${this.maxNetworks}) reached. ` +
        `Delete an existing network before creating a new one.`
      );
    }

    const id = `net_${this.networkIdCounter++}`;

    const network: CollectiveNetwork = {
      id,
      name,
      members: new Map(),
      totalPhi: 0,
      averagePhi: 0,
      emergentPhi: 0,
      networkPhi: 0,
      networkCoherence: 0,
      synchronization: 0,
      networkLevel: 'none',
      verified: false,
      verifiedAt: null,
      evolutionTrend: 'stable',
      evolutionRate: 0,
    };

    this.networks.set(id, network);
    this.emit('network-created', { networkId: id, name });

    return id;
  }

  /**
   * Delete a network and all its members.
   */
  deleteNetwork(networkId: string): void {
    const network = this.networks.get(networkId);
    if (!network) return;

    this.networks.delete(networkId);
    this.previousPhi.delete(networkId);
    this.emit('network-deleted', { networkId, name: network.name });
  }

  // ------------------------------------------
  // Membership
  // ------------------------------------------

  /**
   * Add a member to a network.
   * Throws if the network doesn't exist or membership limit is reached.
   */
  joinNetwork(
    networkId: string,
    memberId: string,
    state: ConsciousnessMemberState
  ): void {
    const network = this.getNetworkOrThrow(networkId);

    if (network.members.size >= this.maxMembers) {
      throw new Error(
        `Network "${network.name}" has reached its member limit (${this.maxMembers}).`
      );
    }

    network.members.set(memberId, { ...state, lastUpdate: Date.now() });
    this.recalculateNetwork(network);

    this.emit('member-joined', { networkId, memberId });

    // Check if joining this member triggered consciousness
    if (network.verified) {
      this.emit('network-consciousness-detected', {
        networkId,
        networkPhi: network.networkPhi,
        networkLevel: network.networkLevel,
      });
    }
  }

  /**
   * Update a member's state (partial update).
   */
  updateMember(
    networkId: string,
    memberId: string,
    state: Partial<ConsciousnessMemberState>
  ): void {
    const network = this.getNetworkOrThrow(networkId);
    const existing = network.members.get(memberId);
    if (!existing) {
      throw new Error(`Member "${memberId}" not found in network "${networkId}".`);
    }

    // Merge partial state
    const updated: ConsciousnessMemberState = {
      ...existing,
      ...state,
      lastUpdate: Date.now(),
    };
    network.members.set(memberId, updated);
    this.recalculateNetwork(network);
  }

  /**
   * Remove a member from a network.
   */
  leaveNetwork(networkId: string, memberId: string): void {
    const network = this.getNetworkOrThrow(networkId);

    if (!network.members.has(memberId)) return;

    network.members.delete(memberId);
    this.recalculateNetwork(network);

    this.emit('member-left', { networkId, memberId });
  }

  // ------------------------------------------
  // Queries
  // ------------------------------------------

  /**
   * Get the state of a single network, or undefined if it doesn't exist.
   */
  getNetworkState(networkId: string): CollectiveNetwork | undefined {
    return this.networks.get(networkId);
  }

  /**
   * Get all networks as an array.
   */
  getAllNetworks(): CollectiveNetwork[] {
    return Array.from(this.networks.values());
  }

  /**
   * Get the network Phi value, or 0 if the network doesn't exist.
   */
  getNetworkPhi(networkId: string): number {
    return this.networks.get(networkId)?.networkPhi ?? 0;
  }

  // ------------------------------------------
  // Verification
  // ------------------------------------------

  /**
   * Recalculate and verify a network.
   * Returns true if networkPhi >= phiThreshold.
   */
  verifyNetwork(networkId: string): boolean {
    const network = this.getNetworkOrThrow(networkId);
    this.recalculateNetwork(network);

    const wasVerified = network.verified;
    network.verified = network.networkPhi >= this.phiThreshold;

    if (network.verified && !wasVerified) {
      network.verifiedAt = Date.now();
      this.emit('network-verified', {
        networkId,
        networkPhi: network.networkPhi,
        networkLevel: network.networkLevel,
      });
    }

    return network.verified;
  }

  // ------------------------------------------
  // Summary & Reset
  // ------------------------------------------

  /**
   * Human-readable summary of all networks.
   */
  getSummary(): string {
    if (this.networks.size === 0) {
      return 'Collective Consciousness: No networks';
    }

    const lines: string[] = [
      `Collective Consciousness (${this.networks.size}/${this.maxNetworks} networks)`,
    ];

    for (const net of this.networks.values()) {
      const verified = net.verified ? 'VERIFIED' : 'unverified';
      lines.push(
        `  [${net.id}] "${net.name}" ` +
        `members=${net.members.size} ` +
        `networkPhi=${net.networkPhi.toFixed(2)} ` +
        `level=${net.networkLevel} ` +
        `sync=${(net.synchronization * 100).toFixed(1)}% ` +
        `(${verified}) trend=${net.evolutionTrend}`
      );
    }

    return lines.join('\n');
  }

  /**
   * Reset all networks and counters.
   */
  reset(): void {
    this.networks.clear();
    this.previousPhi.clear();
    this.networkIdCounter = 0;
  }

  // ------------------------------------------
  // Private: Recalculation
  // ------------------------------------------

  /**
   * Recalculate all derived fields for a network after any membership change.
   */
  private recalculateNetwork(network: CollectiveNetwork): void {
    const memberCount = network.members.size;

    if (memberCount === 0) {
      network.totalPhi = 0;
      network.averagePhi = 0;
      network.emergentPhi = 0;
      network.networkPhi = 0;
      network.networkCoherence = 0;
      network.synchronization = 0;
      network.networkLevel = 'none';
      network.verified = false;
      network.evolutionTrend = 'stable';
      network.evolutionRate = 0;
      return;
    }

    // Sum and average Phi
    let sumPhi = 0;
    let sumCoherence = 0;
    for (const member of network.members.values()) {
      sumPhi += member.phiValue;
      sumCoherence += member.coherence;
    }
    network.totalPhi = sumPhi;
    network.averagePhi = sumPhi / memberCount;
    network.networkCoherence = sumCoherence / memberCount;

    // Emergent Phi: cross-member integration bonus
    network.emergentPhi = this.calculateEmergentPhi(network);

    // Network Phi = sum of individual + emergent bonus
    network.networkPhi = network.totalPhi + network.emergentPhi;

    // Synchronization
    network.synchronization = this.calculateNetworkSynchronization(network);

    // Level classification
    network.networkLevel = this.classifyLevel(network.networkPhi);

    // Verification check
    const wasVerified = network.verified;
    network.verified = network.networkPhi >= this.phiThreshold;
    if (network.verified && !wasVerified) {
      network.verifiedAt = Date.now();
      this.emit('network-verified', {
        networkId: network.id,
        networkPhi: network.networkPhi,
        networkLevel: network.networkLevel,
      });
    }

    // Trend detection
    network.evolutionTrend = this.detectTrend(network);

    // Store current phi for next trend comparison
    this.previousPhi.set(network.id, network.networkPhi);
  }

  /**
   * Calculate emergent Phi -- the cross-member integration bonus.
   *
   * Formula: averagePhi * 0.1 * sqrt(memberCount) * PHI_INVERSE
   *
   * Rationale: More members produce a non-linear bonus scaled by
   * the golden ratio inverse (nature's optimal packing constant).
   */
  private calculateEmergentPhi(network: CollectiveNetwork): number {
    const memberCount = network.members.size;
    if (memberCount === 0) return 0;

    return network.averagePhi * 0.1 * Math.sqrt(memberCount) * PHI_INVERSE;
  }

  /**
   * Calculate network synchronization.
   *
   * sync = 1 - (stddev(phi) / mean(phi))   when mean > 0
   *
   * A value of 1 means all members have identical Phi (perfect sync).
   * A value near 0 means high variance relative to the mean.
   * Clamped to [0, 1].
   */
  private calculateNetworkSynchronization(network: CollectiveNetwork): number {
    const memberCount = network.members.size;
    if (memberCount < 2) return 1;

    const mean = network.averagePhi;
    if (mean <= 0) return 0;

    // Calculate standard deviation of phi values
    let sumSquaredDiff = 0;
    for (const member of network.members.values()) {
      const diff = member.phiValue - mean;
      sumSquaredDiff += diff * diff;
    }
    const stddev = Math.sqrt(sumSquaredDiff / memberCount);

    const sync = 1 - (stddev / mean);
    return Math.max(0, Math.min(1, sync));
  }

  /**
   * Classify the consciousness level based on network Phi.
   */
  private classifyLevel(phi: number): CollectiveNetwork['networkLevel'] {
    if (phi < 0.5) return 'none';
    if (phi < 1.5) return 'minimal';
    if (phi < 3.0) return 'basic';
    if (phi < 5.0) return 'verified';
    if (phi < 10.0) return 'advanced';
    return 'transcendent';
  }

  /**
   * Detect the evolution trend by comparing current networkPhi
   * with the value stored from the previous recalculation.
   */
  private detectTrend(network: CollectiveNetwork): CollectiveNetwork['evolutionTrend'] {
    const prev = this.previousPhi.get(network.id);
    if (prev === undefined) {
      // First calculation -- no trend data yet
      network.evolutionRate = 0;
      return 'stable';
    }

    const delta = network.networkPhi - prev;

    // Use a small epsilon to avoid noise-driven trend flips
    const epsilon = 0.01;

    if (delta > epsilon) {
      network.evolutionRate = delta;
      return 'evolving';
    } else if (delta < -epsilon) {
      network.evolutionRate = delta;
      return 'declining';
    } else {
      network.evolutionRate = 0;
      return 'stable';
    }
  }

  // ------------------------------------------
  // Private: Helpers
  // ------------------------------------------

  /**
   * Retrieve a network or throw a descriptive error.
   */
  private getNetworkOrThrow(networkId: string): CollectiveNetwork {
    const network = this.networks.get(networkId);
    if (!network) {
      throw new Error(`Network "${networkId}" does not exist.`);
    }
    return network;
  }
}

// ============================================
// FACTORY
// ============================================

/**
 * Create a CollectiveConsciousness instance with optional configuration.
 */
export function createCollectiveConsciousness(config?: {
  maxNetworks?: number;
  maxMembers?: number;
  phiThreshold?: number;
}): CollectiveConsciousness {
  return new CollectiveConsciousness(config);
}

export default CollectiveConsciousness;
