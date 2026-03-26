#pragma once
#include "../../core/QuantumAlgorithm.h"

namespace QSE {

// Quantum Fourier Transform (QFT)
// Transforms computational basis states using DFT coefficients.
// The QFT on n qubits maps:
//   |j> -> (1/√N) Σ_k e^(2πijk/N) |k>
class QFT : public QuantumAlgorithm {
public:
    // If inverse=true, applies QFT†
    QFT(int numQubits, bool inverse = false);

    void buildCircuit() override;
    void run() override;
    std::string resultSummary() const override;

    // Apply QFT to specific qubit range [startQubit, startQubit+length)
    // within an existing register
    static void apply(QubitRegister& qreg, int startQubit, int length, bool inverse = false);
    static void applyFull(QubitRegister& qreg, bool inverse = false);

private:
    bool inverse_;
    void applyQFT(QubitRegister& qreg, int start, int n, bool inv);
};

} // namespace QSE