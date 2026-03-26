#pragma once
#include "../../core/QuantumAlgorithm.h"
#include <vector>
#include <random>
#include <string>

namespace QSE {

class QMeans : public QuantumAlgorithm {
public:
    QMeans(int numQubits, int k, unsigned seed = 42);

    void buildCircuit() override {}
    void run() override;
    std::string resultSummary() const override { return description(); }
    std::string description() const;

    void setData(const std::vector<std::vector<double>>& data);
    void setCentroids(const std::vector<std::vector<double>>& centroids);

    std::vector<int>                  getLabels()     const { return m_labels; }
    std::vector<std::vector<double>>  getCentroids()  const { return m_centroids; }
    int                               getIterations() const { return m_itersRun; }

    double swapTestSimilarity(const std::vector<double>& a,
                              const std::vector<double>& b,
                              int shots = 1000);

private:
    int                              m_k;
    std::vector<std::vector<double>> m_data;
    std::vector<std::vector<double>> m_centroids;
    std::vector<int>                 m_labels;
    int                              m_itersRun;
    std::mt19937                     m_rng;

    void initCentroidsRandom();
    void assignLabels();
    void updateCentroids();
    bool converged(const std::vector<std::vector<double>>& oldC) const;
    double quantumDistance(const std::vector<double>& a,
                           const std::vector<double>& b);
    void amplitudeEncode(QubitRegister& reg, const std::vector<double>& v,
                         int startQubit, int numQubits);
};

} // namespace QSE