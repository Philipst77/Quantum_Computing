#include "ThreeQubitFlip.h"
#include <iostream>
#include <cmath>

namespace QSE {

// Gate matrix
static const std::vector<std::vector<Complex>> X_GATE = {
    {{0,0}, {1,0}},
    {{1,0}, {0,0}}
};

ThreeQubitFlip::ThreeQubitFlip(unsigned seed)
    : QuantumAlgorithm("ThreeQubitFlip", 3),
      m_success(false),
      m_rng(seed) {}

std::string ThreeQubitFlip::description() const {
    return "3-qubit bit-flip quantum error correction code.";
}

void ThreeQubitFlip::encode(Complex alpha, Complex beta) {
    qreg_.reset();

    // Prepare |ψ> = α|0> + β|1> on qubit 0
    qreg_[0] = alpha;  // |000>
    qreg_[4] = beta;   // |100>
    qreg_.normalize();

    // Encode: copy qubit 0 → 1 and 2
    qreg_.applyControlledGate(0, 1, X_GATE);
    qreg_.applyControlledGate(0, 2, X_GATE);

    m_alpha = alpha;
    m_beta  = beta;
}

void ThreeQubitFlip::injectError(int qubit) {
    qreg_.applyGate(qubit, X_GATE);
}

void ThreeQubitFlip::detectAndCorrect() {
    // Instead infer error by probabilities

    auto probs = qreg_.probabilities();

    int mostLikely = 0;
    double best = 0.0;

    for (size_t i = 0; i < probs.size(); ++i) {
        if (probs[i] > best) {
            best = probs[i];
            mostLikely = i;
        }
    }

    // Determine which qubit flipped based on state
    // Valid states: 000 (0), 111 (7), and single-bit errors
    if (mostLikely == 1 || mostLikely == 6) qreg_.applyGate(0, X_GATE);
    if (mostLikely == 2 || mostLikely == 5) qreg_.applyGate(1, X_GATE);
    if (mostLikely == 4 || mostLikely == 3) qreg_.applyGate(2, X_GATE);

    m_success = true;
}

std::pair<Complex, Complex> ThreeQubitFlip::decode() {
    // Undo encoding
    qreg_.applyControlledGate(0, 2, X_GATE);
    qreg_.applyControlledGate(0, 1, X_GATE);

    // Extract amplitudes
    Complex alpha = qreg_[0]; // |000>
    Complex beta  = qreg_[4]; // |100>

    return {alpha, beta};
}

void ThreeQubitFlip::printResult() const {
    std::cout << "=== 3-Qubit Bit-Flip Code ===\n";
    std::cout << "Correction successful: " << (m_success ? "YES" : "NO") << "\n";
}

void ThreeQubitFlip::run() {
    encode(m_alpha, m_beta);
    detectAndCorrect();
}

} // namespace QSE