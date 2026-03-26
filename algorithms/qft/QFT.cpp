#include "QFT.h"
#include <sstream>
#include <cmath>

namespace QSE {

QFT::QFT(int numQubits, bool inverse)
    : QuantumAlgorithm("QFT", numQubits), inverse_(inverse) {}

void QFT::applyQFT(QubitRegister& qreg, int start, int n, bool inv) {

    if (!inv) {
        for (int i = start; i < start + n; ++i) {
            for (int j = start; j < i; ++j) {
                double angle = M_PI / (1ULL << (i - j));
                Gates::applyCP(qreg, j, i, angle);
            }
            Gates::applyH(qreg, i);
        }

        for (int i = 0; i < n / 2; ++i)
            Gates::applySWAP(qreg, start + i, start + n - 1 - i);

    } else {

        for (int i = 0; i < n / 2; ++i)
            Gates::applySWAP(qreg, start + i, start + n - 1 - i);

        for (int i = start + n - 1; i >= start; --i) {
            Gates::applyH(qreg, i);
            for (int j = i - 1; j >= start; --j) {
                double angle = -M_PI / (1ULL << (i - j));
                Gates::applyCP(qreg, j, i, angle);
            }
        }
    }
}

void QFT::buildCircuit() {
    qreg_.reset();
    log("QFT circuit built for " + std::to_string(numQubits_) + " qubits");
}

void QFT::run() {
    applyQFT(qreg_, 0, numQubits_, inverse_);
    log("QFT applied (" + std::string(inverse_ ? "inverse" : "forward") + ")");
}

std::string QFT::resultSummary() const {
    std::ostringstream oss;
    oss << "=== QFT Result ===\n"
        << "  Qubits:    " << numQubits_ << "\n"
        << "  Direction: " << (inverse_ ? "Inverse" : "Forward") << "\n"
        << "  State:     " << qreg_.toString() << "\n";
    return oss.str();
}

void QFT::apply(QubitRegister& qreg, int startQubit, int length, bool inverse) {
    QFT helper(length, inverse);
    helper.applyQFT(qreg, startQubit, length, inverse);
}

void QFT::applyFull(QubitRegister& qreg, bool inverse) {
    QFT helper(qreg.numQubits(), inverse);
    helper.applyQFT(qreg, 0, qreg.numQubits(), inverse);
}

} // namespace QSE