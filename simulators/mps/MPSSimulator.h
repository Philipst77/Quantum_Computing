#pragma once
#include "../../core/QubitRegister.h"
#include <vector>
#include <string>
#include <random>

namespace QSE {

class MPSSimulator {
public:
    explicit MPSSimulator(int numQubits, int maxBondDim = 64);

    // Single-qubit gates
    void applyH(int qubit);
    void applyX(int qubit);
    void applyY(int qubit);
    void applyZ(int qubit);
    void applyS(int qubit);
    void applyT(int qubit);
    void applyRz(int qubit, double theta);
    void applyRy(int qubit, double theta);

    // Two-qubit gates (adjacent qubits only for CNOT)
    void applyCNOT(int control, int target);
    void applyTwoSiteCNOT(int control, int target); // arbitrary distance

    // Measurement
    int              measure(int qubit, std::mt19937& rng);
    std::vector<int> measureAll(std::mt19937& rng);

    // State access
    std::vector<Complex> toStateVector() const;
    double               probability(int qubit, int outcome) const;

    // Utilities
    int              numQubits()    const { return m_n; }
    int              maxBondDim()   const { return m_chi; }
    std::vector<int> bondDimensions() const;
    void             reset();

private:
    int m_n;
    int m_chi;

    // Tensor shape: [physical_dim=2][chiL][chiR]
    using Tensor3 = std::vector<std::vector<std::vector<Complex>>>;
    std::vector<Tensor3>               m_tensors;
    std::vector<std::vector<double>>   m_sv; // singular values per bond

    Tensor3 makeTensor(int d, int chiL, int chiR) const;
    void    applySingleSiteGate(int q, const std::vector<std::vector<Complex>>& gate);
    void    applyTwoSiteGate(int q, const std::vector<std::vector<Complex>>& gate);
    void    svdAndTruncate(int bond);

    struct SVDResult {
        std::vector<std::vector<Complex>> U, Vt;
        std::vector<double>               S;
    };
    SVDResult thinSVD(const std::vector<std::vector<Complex>>& M, int maxRank) const;
};

} // namespace QSE