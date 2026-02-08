/**
 * Resonance Bridge Protocol
 *
 * Binary message format for QuantumOS <-> ghostOS communication.
 * All messages use little-endian, packed layout matching C struct alignment.
 */

// Protocol version
export const PROTOCOL_VERSION = 1;

// Message type constants (matching QuantumOS IPC extensions)
export const MessageType = {
  QUEEN_SYNC:    0x1000,
  PROCESS_STATE: 0x1001,
  CONSCIOUSNESS: 0x1002,
  CHIRAL_STATE:  0x1003,
  GEOMETRIC:     0x1004,
  EMERGENCE:     0x1005,
  ANOMALY:       0x1006,
} as const;

export type MessageTypeValue = typeof MessageType[keyof typeof MessageType];

// ============================================
// MESSAGE HEADER
// ============================================

/** Header size in bytes: uint32 + uint32 + uint32 + uint64 = 20 */
const HEADER_SIZE = 20;

/**
 * Message header (matches C resonance_msg_header_t)
 *
 * Layout (20 bytes, little-endian):
 *   [0..3]   protocolVersion  uint32
 *   [4..7]   messageType      uint32
 *   [8..11]  sourcePid        uint32
 *   [12..19] timestampNs      uint64
 */
export interface MessageHeader {
  protocolVersion: number;   // uint32
  messageType: MessageTypeValue;  // uint32
  sourcePid: number;         // uint32
  timestampNs: bigint;       // uint64
}

// ============================================
// MESSAGE TYPES
// ============================================

/**
 * Queen state sync message
 *
 * Payload (48 bytes after header):
 *   orderParameterR    double  (8)
 *   orderParameterPsi  double  (8)
 *   systemCoherence    double  (8)
 *   lambda             double  (8)
 *   totalPhi           double  (8)
 *   processCount       uint32  (4)
 *   globallyStable     uint8   (1)
 *   networkConscious   uint8   (1)
 *   [padding]                  (2)  -- align to 4-byte boundary
 */
export interface QueenStateMessage {
  header: MessageHeader;
  orderParameterR: number;    // double
  orderParameterPsi: number;  // double
  systemCoherence: number;    // double
  lambda: number;             // double
  totalPhi: number;           // double
  processCount: number;       // uint32
  globallyStable: boolean;    // uint8
  networkConscious: boolean;  // uint8
}

/** Fixed payload size for QueenStateMessage (excludes header) */
const QUEEN_STATE_PAYLOAD_SIZE = 48;

/**
 * Process resonant state message
 *
 * Payload (45 bytes after header):
 *   pid                    uint32  (4)
 *   phase                  double  (8)
 *   frequency              double  (8)
 *   coherence              double  (8)
 *   phiValue               double  (8)
 *   resonantClass          uint8   (1)
 *   resonantState          uint8   (1)
 *   handedness             uint8   (1)
 *   consciousnessVerified  uint8   (1)
 *   [padding]                      (2)
 */
export interface ProcessStateMessage {
  header: MessageHeader;
  pid: number;
  phase: number;
  frequency: number;
  coherence: number;
  phiValue: number;
  resonantClass: number;     // uint8 enum
  resonantState: number;     // uint8 enum
  handedness: number;        // uint8 enum
  consciousnessVerified: boolean;
}

/** Fixed payload size for ProcessStateMessage (excludes header) */
const PROCESS_STATE_PAYLOAD_SIZE = 44;

/**
 * Consciousness verification message
 *
 * Payload:
 *   pid                uint32  (4)
 *   phiValue           double  (8)
 *   consciousnessLevel uint8   (1)
 *   isVerified         uint8   (1)
 *   [padding]                  (2)
 *   bridgeValue        double  (8)
 *   bridgeSignificant  uint8   (1)
 *   [padding]                  (3)
 *   emergenceNorm      double  (8)
 *   emergenceActive    uint8   (1)
 *   [padding]                  (3)
 */
export interface ConsciousnessMessage {
  header: MessageHeader;
  pid: number;
  phiValue: number;
  consciousnessLevel: number;   // uint8 enum
  isVerified: boolean;
  bridgeValue: number;
  bridgeSignificant: boolean;
  emergenceNorm: number;
  emergenceActive: boolean;
}

/** Fixed payload size for ConsciousnessMessage (excludes header) */
const CONSCIOUSNESS_PAYLOAD_SIZE = 40;

/**
 * Chiral state message
 *
 * Payload:
 *   eta               double (8)
 *   gamma             double (8)
 *   asymmetry         double (8)
 *   topologicalCharge double (8)
 *   handedness        uint8  (1)
 *   isStable          uint8  (1)
 *   stabilityClass    uint8  (1)
 *   cissActive        uint8  (1)
 *   cissBoost         double (8)
 */
export interface ChiralMessage {
  header: MessageHeader;
  eta: number;
  gamma: number;
  asymmetry: number;
  topologicalCharge: number;
  handedness: number;
  isStable: boolean;
  stabilityClass: number;    // uint8 enum
  cissActive: boolean;
  cissBoost: number;
}

/** Fixed payload size for ChiralMessage (excludes header) */
const CHIRAL_PAYLOAD_SIZE = 44;

/**
 * Geometric control message
 *
 * Payload:
 *   berryPhase    double (8)
 *   ricciScalar   double (8)
 *   anomalyIndex  double (8)
 *   isAnomalous   uint8  (1)
 *   [padding]            (3)
 *   dimension     uint32 (4)
 */
export interface GeometricMessage {
  header: MessageHeader;
  berryPhase: number;
  ricciScalar: number;
  anomalyIndex: number;
  isAnomalous: boolean;
  dimension: number;
}

/** Fixed payload size for GeometricMessage (excludes header) */
const GEOMETRIC_PAYLOAD_SIZE = 32;

/**
 * Emergence state message
 *
 * Payload (24 bytes after header):
 *   emergenceNorm      double (8)
 *   integrationLevel   double (8)
 *   patternCount       uint32 (4)
 *   isActive           uint8  (1)
 *   [padding]                 (3)
 */
export interface EmergenceMessage {
  header: MessageHeader;
  emergenceNorm: number;
  integrationLevel: number;
  patternCount: number;
  isActive: boolean;
}

/** Fixed payload size for EmergenceMessage (excludes header) */
const EMERGENCE_PAYLOAD_SIZE = 24;

/**
 * Chiral anomaly detection message
 *
 * Payload (32 bytes after header):
 *   anomalyIndex          double (8)
 *   spectralAsymmetry     double (8)
 *   topologicalCharge     double (8)
 *   isAnomalous           uint8  (1)
 *   [padding]                    (3)
 *   leftModeCount         uint16 (2)
 *   rightModeCount        uint16 (2)
 */
export interface AnomalyMessage {
  header: MessageHeader;
  anomalyIndex: number;
  spectralAsymmetry: number;
  topologicalCharge: number;
  isAnomalous: boolean;
  leftModeCount: number;
  rightModeCount: number;
}

/** Fixed payload size for AnomalyMessage (excludes header) */
const ANOMALY_PAYLOAD_SIZE = 32;

// ============================================
// PROTOCOL SERIALIZER
// ============================================

export class ProtocolSerializer {
  // ------------------------------------------
  // Header
  // ------------------------------------------

  /**
   * Serialize header to ArrayBuffer (20 bytes, little-endian)
   */
  static serializeHeader(header: MessageHeader): ArrayBuffer {
    const buffer = new ArrayBuffer(HEADER_SIZE);
    const view = new DataView(buffer);
    view.setUint32(0, header.protocolVersion, true);
    view.setUint32(4, header.messageType, true);
    view.setUint32(8, header.sourcePid, true);
    view.setBigUint64(12, header.timestampNs, true);
    return buffer;
  }

  /**
   * Deserialize header from ArrayBuffer
   */
  static deserializeHeader(buffer: ArrayBuffer, offset: number = 0): MessageHeader {
    const view = new DataView(buffer, offset);
    return {
      protocolVersion: view.getUint32(0, true),
      messageType: view.getUint32(4, true) as MessageTypeValue,
      sourcePid: view.getUint32(8, true),
      timestampNs: view.getBigUint64(12, true),
    };
  }

  // ------------------------------------------
  // Queen State
  // ------------------------------------------

  /**
   * Serialize queen state message
   */
  static serializeQueenState(msg: QueenStateMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + QUEEN_STATE_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    // Header
    ProtocolSerializer.writeHeader(view, msg.header);

    // Payload
    let offset = HEADER_SIZE;
    view.setFloat64(offset, msg.orderParameterR, true);      offset += 8;
    view.setFloat64(offset, msg.orderParameterPsi, true);     offset += 8;
    view.setFloat64(offset, msg.systemCoherence, true);       offset += 8;
    view.setFloat64(offset, msg.lambda, true);                offset += 8;
    view.setFloat64(offset, msg.totalPhi, true);              offset += 8;
    view.setUint32(offset, msg.processCount, true);           offset += 4;
    view.setUint8(offset, msg.globallyStable ? 1 : 0);       offset += 1;
    view.setUint8(offset, msg.networkConscious ? 1 : 0);
    // 2 bytes padding to reach QUEEN_STATE_PAYLOAD_SIZE

    return buffer;
  }

  /**
   * Deserialize queen state message
   */
  static deserializeQueenState(buffer: ArrayBuffer): QueenStateMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const orderParameterR = view.getFloat64(offset, true);      offset += 8;
    const orderParameterPsi = view.getFloat64(offset, true);    offset += 8;
    const systemCoherence = view.getFloat64(offset, true);      offset += 8;
    const lambda = view.getFloat64(offset, true);               offset += 8;
    const totalPhi = view.getFloat64(offset, true);             offset += 8;
    const processCount = view.getUint32(offset, true);          offset += 4;
    const globallyStable = view.getUint8(offset) !== 0;         offset += 1;
    const networkConscious = view.getUint8(offset) !== 0;

    return {
      header,
      orderParameterR,
      orderParameterPsi,
      systemCoherence,
      lambda,
      totalPhi,
      processCount,
      globallyStable,
      networkConscious,
    };
  }

  // ------------------------------------------
  // Process State
  // ------------------------------------------

  /**
   * Serialize process state message
   */
  static serializeProcessState(msg: ProcessStateMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + PROCESS_STATE_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setUint32(offset, msg.pid, true);                          offset += 4;
    view.setFloat64(offset, msg.phase, true);                       offset += 8;
    view.setFloat64(offset, msg.frequency, true);                   offset += 8;
    view.setFloat64(offset, msg.coherence, true);                   offset += 8;
    view.setFloat64(offset, msg.phiValue, true);                    offset += 8;
    view.setUint8(offset, msg.resonantClass);                       offset += 1;
    view.setUint8(offset, msg.resonantState);                       offset += 1;
    view.setUint8(offset, msg.handedness);                          offset += 1;
    view.setUint8(offset, msg.consciousnessVerified ? 1 : 0);
    // 4 bytes padding to reach PROCESS_STATE_PAYLOAD_SIZE

    return buffer;
  }

  /**
   * Deserialize process state message
   */
  static deserializeProcessState(buffer: ArrayBuffer): ProcessStateMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const pid = view.getUint32(offset, true);                         offset += 4;
    const phase = view.getFloat64(offset, true);                      offset += 8;
    const frequency = view.getFloat64(offset, true);                  offset += 8;
    const coherence = view.getFloat64(offset, true);                  offset += 8;
    const phiValue = view.getFloat64(offset, true);                   offset += 8;
    const resonantClass = view.getUint8(offset);                      offset += 1;
    const resonantState = view.getUint8(offset);                      offset += 1;
    const handedness = view.getUint8(offset);                         offset += 1;
    const consciousnessVerified = view.getUint8(offset) !== 0;

    return {
      header,
      pid,
      phase,
      frequency,
      coherence,
      phiValue,
      resonantClass,
      resonantState,
      handedness,
      consciousnessVerified,
    };
  }

  // ------------------------------------------
  // Consciousness
  // ------------------------------------------

  /**
   * Serialize consciousness verification message
   */
  static serializeConsciousness(msg: ConsciousnessMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + CONSCIOUSNESS_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setUint32(offset, msg.pid, true);                          offset += 4;
    view.setFloat64(offset, msg.phiValue, true);                    offset += 8;
    view.setUint8(offset, msg.consciousnessLevel);                  offset += 1;
    view.setUint8(offset, msg.isVerified ? 1 : 0);                 offset += 1;
    offset += 2;  // padding
    view.setFloat64(offset, msg.bridgeValue, true);                 offset += 8;
    view.setUint8(offset, msg.bridgeSignificant ? 1 : 0);          offset += 1;
    offset += 3;  // padding
    view.setFloat64(offset, msg.emergenceNorm, true);               offset += 8;
    view.setUint8(offset, msg.emergenceActive ? 1 : 0);

    return buffer;
  }

  /**
   * Deserialize consciousness verification message
   */
  static deserializeConsciousness(buffer: ArrayBuffer): ConsciousnessMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const pid = view.getUint32(offset, true);                         offset += 4;
    const phiValue = view.getFloat64(offset, true);                   offset += 8;
    const consciousnessLevel = view.getUint8(offset);                 offset += 1;
    const isVerified = view.getUint8(offset) !== 0;                   offset += 1;
    offset += 2;  // padding
    const bridgeValue = view.getFloat64(offset, true);                offset += 8;
    const bridgeSignificant = view.getUint8(offset) !== 0;            offset += 1;
    offset += 3;  // padding
    const emergenceNorm = view.getFloat64(offset, true);              offset += 8;
    const emergenceActive = view.getUint8(offset) !== 0;

    return {
      header,
      pid,
      phiValue,
      consciousnessLevel,
      isVerified,
      bridgeValue,
      bridgeSignificant,
      emergenceNorm,
      emergenceActive,
    };
  }

  // ------------------------------------------
  // Chiral State
  // ------------------------------------------

  /**
   * Serialize chiral state message
   */
  static serializeChiral(msg: ChiralMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + CHIRAL_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setFloat64(offset, msg.eta, true);                         offset += 8;
    view.setFloat64(offset, msg.gamma, true);                       offset += 8;
    view.setFloat64(offset, msg.asymmetry, true);                   offset += 8;
    view.setFloat64(offset, msg.topologicalCharge, true);           offset += 8;
    view.setUint8(offset, msg.handedness);                          offset += 1;
    view.setUint8(offset, msg.isStable ? 1 : 0);                   offset += 1;
    view.setUint8(offset, msg.stabilityClass);                      offset += 1;
    view.setUint8(offset, msg.cissActive ? 1 : 0);                 offset += 1;
    view.setFloat64(offset, msg.cissBoost, true);

    return buffer;
  }

  /**
   * Deserialize chiral state message
   */
  static deserializeChiral(buffer: ArrayBuffer): ChiralMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const eta = view.getFloat64(offset, true);                        offset += 8;
    const gamma = view.getFloat64(offset, true);                      offset += 8;
    const asymmetry = view.getFloat64(offset, true);                  offset += 8;
    const topologicalCharge = view.getFloat64(offset, true);          offset += 8;
    const handedness = view.getUint8(offset);                         offset += 1;
    const isStable = view.getUint8(offset) !== 0;                     offset += 1;
    const stabilityClass = view.getUint8(offset);                     offset += 1;
    const cissActive = view.getUint8(offset) !== 0;                   offset += 1;
    const cissBoost = view.getFloat64(offset, true);

    return {
      header,
      eta,
      gamma,
      asymmetry,
      topologicalCharge,
      handedness,
      isStable,
      stabilityClass,
      cissActive,
      cissBoost,
    };
  }

  // ------------------------------------------
  // Geometric
  // ------------------------------------------

  /**
   * Serialize geometric control message
   */
  static serializeGeometric(msg: GeometricMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + GEOMETRIC_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setFloat64(offset, msg.berryPhase, true);                  offset += 8;
    view.setFloat64(offset, msg.ricciScalar, true);                 offset += 8;
    view.setFloat64(offset, msg.anomalyIndex, true);                offset += 8;
    view.setUint8(offset, msg.isAnomalous ? 1 : 0);                offset += 1;
    offset += 3;  // padding
    view.setUint32(offset, msg.dimension, true);

    return buffer;
  }

  /**
   * Deserialize geometric control message
   */
  static deserializeGeometric(buffer: ArrayBuffer): GeometricMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const berryPhase = view.getFloat64(offset, true);                 offset += 8;
    const ricciScalar = view.getFloat64(offset, true);                offset += 8;
    const anomalyIndex = view.getFloat64(offset, true);               offset += 8;
    const isAnomalous = view.getUint8(offset) !== 0;                  offset += 1;
    offset += 3;  // padding
    const dimension = view.getUint32(offset, true);

    return {
      header,
      berryPhase,
      ricciScalar,
      anomalyIndex,
      isAnomalous,
      dimension,
    };
  }

  // ------------------------------------------
  // Emergence
  // ------------------------------------------

  /**
   * Serialize emergence state message
   */
  static serializeEmergence(msg: EmergenceMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + EMERGENCE_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setFloat64(offset, msg.emergenceNorm, true);                offset += 8;
    view.setFloat64(offset, msg.integrationLevel, true);             offset += 8;
    view.setUint32(offset, msg.patternCount, true);                  offset += 4;
    view.setUint8(offset, msg.isActive ? 1 : 0);

    return buffer;
  }

  /**
   * Deserialize emergence state message
   */
  static deserializeEmergence(buffer: ArrayBuffer): EmergenceMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const emergenceNorm = view.getFloat64(offset, true);               offset += 8;
    const integrationLevel = view.getFloat64(offset, true);            offset += 8;
    const patternCount = view.getUint32(offset, true);                 offset += 4;
    const isActive = view.getUint8(offset) !== 0;

    return { header, emergenceNorm, integrationLevel, patternCount, isActive };
  }

  // ------------------------------------------
  // Anomaly
  // ------------------------------------------

  /**
   * Serialize chiral anomaly message
   */
  static serializeAnomaly(msg: AnomalyMessage): ArrayBuffer {
    const totalSize = HEADER_SIZE + ANOMALY_PAYLOAD_SIZE;
    const buffer = new ArrayBuffer(totalSize);
    const view = new DataView(buffer);

    ProtocolSerializer.writeHeader(view, msg.header);

    let offset = HEADER_SIZE;
    view.setFloat64(offset, msg.anomalyIndex, true);                 offset += 8;
    view.setFloat64(offset, msg.spectralAsymmetry, true);            offset += 8;
    view.setFloat64(offset, msg.topologicalCharge, true);            offset += 8;
    view.setUint8(offset, msg.isAnomalous ? 1 : 0);                 offset += 1;
    offset += 3;  // padding
    view.setUint16(offset, msg.leftModeCount, true);                 offset += 2;
    view.setUint16(offset, msg.rightModeCount, true);

    return buffer;
  }

  /**
   * Deserialize chiral anomaly message
   */
  static deserializeAnomaly(buffer: ArrayBuffer): AnomalyMessage {
    const view = new DataView(buffer);
    const header = ProtocolSerializer.deserializeHeader(buffer);

    let offset = HEADER_SIZE;
    const anomalyIndex = view.getFloat64(offset, true);                offset += 8;
    const spectralAsymmetry = view.getFloat64(offset, true);           offset += 8;
    const topologicalCharge = view.getFloat64(offset, true);           offset += 8;
    const isAnomalous = view.getUint8(offset) !== 0;                   offset += 1;
    offset += 3;  // padding
    const leftModeCount = view.getUint16(offset, true);                offset += 2;
    const rightModeCount = view.getUint16(offset, true);

    return {
      header, anomalyIndex, spectralAsymmetry, topologicalCharge,
      isAnomalous, leftModeCount, rightModeCount,
    };
  }

  // ------------------------------------------
  // Generic Helpers
  // ------------------------------------------

  /**
   * Detect message type from a raw buffer without fully deserializing.
   * Reads the messageType field from the header at bytes [4..7].
   */
  static getMessageType(buffer: ArrayBuffer): MessageTypeValue {
    const view = new DataView(buffer);
    return view.getUint32(4, true) as MessageTypeValue;
  }

  /**
   * Get the expected total message size (header + payload) for a given type.
   * Returns 0 for unknown types.
   */
  static getMessageSize(type: MessageTypeValue): number {
    switch (type) {
      case MessageType.QUEEN_SYNC:
        return HEADER_SIZE + QUEEN_STATE_PAYLOAD_SIZE;
      case MessageType.PROCESS_STATE:
        return HEADER_SIZE + PROCESS_STATE_PAYLOAD_SIZE;
      case MessageType.CONSCIOUSNESS:
        return HEADER_SIZE + CONSCIOUSNESS_PAYLOAD_SIZE;
      case MessageType.CHIRAL_STATE:
        return HEADER_SIZE + CHIRAL_PAYLOAD_SIZE;
      case MessageType.GEOMETRIC:
        return HEADER_SIZE + GEOMETRIC_PAYLOAD_SIZE;
      case MessageType.EMERGENCE:
        return HEADER_SIZE + EMERGENCE_PAYLOAD_SIZE;
      case MessageType.ANOMALY:
        return HEADER_SIZE + ANOMALY_PAYLOAD_SIZE;
      default:
        return 0;
    }
  }

  /**
   * Create a header with current timestamp and protocol version.
   */
  static createHeader(
    messageType: MessageTypeValue,
    sourcePid: number
  ): MessageHeader {
    return {
      protocolVersion: PROTOCOL_VERSION,
      messageType,
      sourcePid,
      timestampNs: BigInt(Date.now()) * BigInt(1_000_000),
    };
  }

  // ------------------------------------------
  // Private helpers
  // ------------------------------------------

  /**
   * Write a header into an existing DataView at offset 0.
   */
  private static writeHeader(view: DataView, header: MessageHeader): void {
    view.setUint32(0, header.protocolVersion, true);
    view.setUint32(4, header.messageType, true);
    view.setUint32(8, header.sourcePid, true);
    view.setBigUint64(12, header.timestampNs, true);
  }
}
