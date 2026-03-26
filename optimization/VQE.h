#pragma once
#include "../core/QuantumAlgorithm.h"
#include <vector>
#include <string>
#include <random>

namespace QSE {

class VQE : public QuantumAlgorithm {
public:
    using PauliTerm = std::pair<double, std::string>;

    VQE(int numQubits, int layers, unsigned seed = 42);

    void buildCircuit() override {}
    void run() override;
    std::string resultSummary() const override { return description(); }
    std::string description() const;

    void setHamiltonian(const std::vector<PauliTerm>& H);
    void setH2Hamiltonian();
    void setParameters(const std::vector<double>& params);

    std::vector<double> getParameters()      const { return m_params; }
    double              getGroundEnergy()    const { return m_groundEnergy; }
    int                 getOptimIterations() const { return m_optIters; }

    double energyExpectation();
    void   printResult() const;

private:
    int                    m_layers;
    std::vector<double>    m_params;
    std::vector<PauliTerm> m_hamiltonian;
    double                 m_groundEnergy;
    int                    m_optIters;
    std::mt19937           m_rng;

    // Rotation gate matrices
    static std::vector<std::vector<Complex>> rxMatrix(double theta);
    static std::vector<std::vector<Complex>> ryMatrix(double theta);
    static std::vector<std::vector<Complex>> rzMatrix(double theta);

    void   applyAnsatz();
    double measurePauliString(const std::string& pauliStr);
    void   gradientDescent(int maxIter, double lr);
    double paramShiftGradient(int paramIdx);
};

} // namespace QSE