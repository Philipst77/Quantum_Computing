#pragma once
#include "../../core/QuantumAlgorithm.h"

namespace QSE {

// Quantum Teleportation Protocol
// Transmits an arbitrary qubit state from Alice to Bob using:
//   - 1 shared EPR pair (Bell state)
//   - 2 classical bits
//
// Qubit layout: [0] = Alice's data qubit, [1] = Alice's EPR half, [2] = Bob's EPR half
class Teleportation : public QuantumAlgorithm {
public:
    // stateToTeleport: [alpha, beta] amplitudes for |0> and |1>
    Teleportation(Complex alpha, Complex beta);

    // Teleport an arbitrary angle (Ry rotation from |0>)
    static Teleportation fromAngle(double theta);

    void buildCircuit() override;
    void run() override;
    std::string resultSummary() const override;

    // After running: retrieve Bob's qubit state fidelity vs original
    double fidelity() const { return fidelity_; }
    bool success() const { return fidelity_ > 0.999; }

private:
    Complex alpha_, beta_;  // Original state amplitudes
    int m1_, m2_;           // Classical measurement results
    double fidelity_;

    void prepareDataQubit();
    void createBellPair();
    void aliceBellMeasurement();
    void bobCorrection();
};

} // namespace QSE