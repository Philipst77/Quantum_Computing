#pragma once

#include "../../core/QuantumAlgorithm.h"
#include <string>
#include <random>

namespace QSE {

class BB84 : public QuantumAlgorithm {
public:
    BB84(int keyLength, unsigned seed = std::random_device{}());

    void buildCircuit() override {}
    void run() override;
    std::string resultSummary() const override;

    void setEavesdropping(bool enable, double interceptProb = 0.5);
    void printResult() const;

    bool isKeySecure()              const { return m_errorRate <= 0.11; }
    const std::string& getSiftedKey() const { return m_siftedKey; }
    double getErrorRate()           const { return m_errorRate; }

private:
    int         m_keyLength;
    bool        m_eavesdrop;
    double      m_interceptProb;
    double      m_errorRate;

    std::mt19937 m_rng;

    std::string m_aliceBits,  m_aliceBases;
    std::string m_bobBits,    m_bobBases;
    std::string m_aliceKey,   m_bobKey;
    std::string m_siftedKey;

    char randomBit();
    char randomBasis();
    void alicePrepare();
    void transmitAndMeasure();
    void siftKeys();
    void estimateError();
};

} // namespace QSE