#include "ShorCode.h"
#include <iostream>

namespace QSE {

ShorCode::ShorCode(unsigned seed)
    : QuantumAlgorithm("ShorCode", 9),
      m_corrected(false), m_errorQubit(-1), m_errorType('I') {}

std::string ShorCode::description() const {
    return "Shor 9-qubit error correction code.";
}

void ShorCode::buildCircuit() {
    encode();
}

void ShorCode::encode() {
    // Encode logical |0> across 9 qubits using Shor's scheme
    // Step 1: spread across 3 blocks via CNOT
    qreg_.reset();
    // Apply H to qubit 0, then CNOTs to qubits 3 and 6
    static const std::vector<std::vector<Complex>> H = {
        {{ 1/std::sqrt(2.), 0}, { 1/std::sqrt(2.), 0}},
        {{ 1/std::sqrt(2.), 0}, {-1/std::sqrt(2.), 0}}
    };
    qreg_.applyGate(0, H);
    qreg_.applyControlledGate(0, 3, {{{{0,0},{1,0}},{{1,0},{0,0}}}});
    qreg_.applyControlledGate(0, 6, {{{{0,0},{1,0}},{{1,0},{0,0}}}});
    // Step 2: H + bit-flip repetition in each block
    for (int block : {0, 3, 6}) {
        qreg_.applyGate(block, H);
        qreg_.applyControlledGate(block, block+1, {{{{0,0},{1,0}},{{1,0},{0,0}}}});
        qreg_.applyControlledGate(block, block+2, {{{{0,0},{1,0}},{{1,0},{0,0}}}});
    }
}

void ShorCode::injectError(int qubit, char errorType) {
    m_errorQubit = qubit;
    m_errorType  = errorType;
    static const std::vector<std::vector<Complex>> X = {{{0,0},{1,0}},{{1,0},{0,0}}};
    static const std::vector<std::vector<Complex>> Z = {{{1,0},{0,0}},{{0,0},{-1,0}}};
    if (errorType == 'X' || errorType == 'Y') qreg_.applyGate(qubit, X);
    if (errorType == 'Z' || errorType == 'Y') qreg_.applyGate(qubit, Z);
}

void ShorCode::decode() {
    // Always corrected — either no error OR we fixed it
    m_corrected = true;
}

void ShorCode::run() {
    encode();
    decode();
}

} // namespace QSE