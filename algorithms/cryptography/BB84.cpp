#include "BB84.h"
#include <iostream>
#include <cmath>
#include <iomanip>

namespace QSE {

// Gate matrices
static const std::vector<std::vector<Complex>> X_GATE = {
    {{0,0}, {1,0}},
    {{1,0}, {0,0}}
};
static const std::vector<std::vector<Complex>> H_GATE = {
    {{ 1.0/std::sqrt(2.0), 0}, { 1.0/std::sqrt(2.0), 0}},
    {{ 1.0/std::sqrt(2.0), 0}, {-1.0/std::sqrt(2.0), 0}}
};

BB84::BB84(int keyLength, unsigned seed)
    : QuantumAlgorithm("BB84", 1),   
      m_keyLength(keyLength),
      m_eavesdrop(false),
      m_interceptProb(0.5),
      m_errorRate(0.0),
      m_rng(seed)
{}

std::string BB84::resultSummary() const {
    return "BB84 QKD: sifted key length=" + std::to_string(m_aliceKey.size()) +
           ", error rate=" + std::to_string(m_errorRate * 100) + "%" +
           ", secure=" + (isKeySecure() ? "YES" : "NO");
}

void BB84::setEavesdropping(bool enable, double interceptProb) {
    m_eavesdrop     = enable;
    m_interceptProb = interceptProb;
}

char BB84::randomBit() {
    std::uniform_int_distribution<int> d(0, 1);
    return '0' + d(m_rng);
}

char BB84::randomBasis() {
    std::uniform_int_distribution<int> d(0, 1);
    return d(m_rng) ? 'X' : 'Z';
}

void BB84::alicePrepare() {
    m_aliceBits.clear();
    m_aliceBases.clear();
    for (int i = 0; i < m_keyLength; ++i) {
        m_aliceBits  += randomBit();
        m_aliceBases += randomBasis();
    }
}

void BB84::transmitAndMeasure() {
    m_bobBases.clear();
    m_bobBits.clear();

    std::uniform_real_distribution<double> prob(0.0, 1.0);

    for (int i = 0; i < m_keyLength; ++i) {
        qreg_.reset();
        if (m_aliceBits[i]  == '1') qreg_.applyGate(0, X_GATE);
        if (m_aliceBases[i] == 'X') qreg_.applyGate(0, H_GATE);

        if (m_eavesdrop && prob(m_rng) < m_interceptProb) {
            char eveBasis = randomBasis();
            if (eveBasis == 'X') qreg_.applyGate(0, H_GATE);
            qreg_.measure(0);
            if (eveBasis == 'X') qreg_.applyGate(0, H_GATE);
        }

        char bobBasis = randomBasis();
        m_bobBases += bobBasis;
        if (bobBasis == 'X') qreg_.applyGate(0, H_GATE);
        int result = qreg_.measure(0);
        m_bobBits += ('0' + result);
    }
}

void BB84::siftKeys() {
    m_aliceKey.clear();
    m_bobKey.clear();
    for (int i = 0; i < m_keyLength; ++i) {
        if (m_aliceBases[i] == m_bobBases[i]) {
            m_aliceKey += m_aliceBits[i];
            m_bobKey   += m_bobBits[i];
        }
    }
    m_siftedKey = m_aliceKey;
}

void BB84::estimateError() {
    if (m_aliceKey.empty()) { m_errorRate = 0.0; return; }
    int errors = 0;
    for (size_t i = 0; i < m_aliceKey.size(); ++i)
        if (m_aliceKey[i] != m_bobKey[i]) ++errors;
    m_errorRate = static_cast<double>(errors) / m_aliceKey.size();
}

void BB84::run() {
    alicePrepare();
    transmitAndMeasure();
    siftKeys();
    estimateError();
}

void BB84::printResult() const {
    std::cout << "=== BB84 QKD Protocol ===\n";
    std::cout << "Raw bits sent    : " << m_keyLength        << "\n";
    std::cout << "Sifted key length: " << m_aliceKey.size()  << "\n";
    std::cout << "Alice sifted key : " << m_aliceKey         << "\n";
    std::cout << "Bob   sifted key : " << m_bobKey           << "\n";
    std::cout << "Error rate       : " << m_errorRate * 100  << "%\n";
    std::cout << "Eavesdropping    : " << (m_eavesdrop ? "YES" : "NO") << "\n";
    std::cout << "Key secure       : " << (isKeySecure() ? "YES" : "NO (>11% QBER)") << "\n";
}

} // namespace QSE