#pragma once
#include "../core/QuantumAlgorithm.h"
#include <string>
#include <vector>

namespace QSE {

// Shor's 9-qubit quantum error correction code
class ShorCode : public QuantumAlgorithm {
public:
    explicit ShorCode(unsigned seed = 42);

    void buildCircuit() override;
    void run() override;
    std::string resultSummary() const override { return description(); }
    std::string description() const;

    void encode();
    void injectError(int qubit, char errorType); // errorType: 'X', 'Z', 'Y'
    void decode();
    bool corrected() const { return m_corrected; }

private:
    bool m_corrected;
    int  m_errorQubit;
    char m_errorType;
};

} // namespace QSE