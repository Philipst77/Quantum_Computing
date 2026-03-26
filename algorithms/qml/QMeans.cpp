#include "QMeans.h"
#include "../../core/Utils.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <algorithm>

namespace QSE {

QMeans::QMeans(int numQubits, int k, unsigned seed)
    : QuantumAlgorithm("QMeans", numQubits), m_k(k), m_itersRun(0), m_rng(seed) {}

std::string QMeans::description() const {
    return "Q-Means: quantum-enhanced k-means clustering with k=" +
           std::to_string(m_k) + " on " + std::to_string(m_data.size()) + " data points.";
}

void QMeans::setData(const std::vector<std::vector<double>>& data) {
    m_data = data;
    m_labels.assign(data.size(), 0);
}

void QMeans::setCentroids(const std::vector<std::vector<double>>& centroids) {
    m_centroids = centroids;
}

void QMeans::amplitudeEncode(QubitRegister& reg, const std::vector<double>& v,
                              int startQubit, int numQubits) {
    size_t dim = 1ULL << numQubits;
    std::vector<Complex> amps(dim, Complex(0,0));
    double norm = 0;
    for (size_t i = 0; i < std::min(v.size(), dim); ++i) norm += v[i] * v[i];
    norm = std::sqrt(norm);
    if (norm < 1e-12) norm = 1.0;
    for (size_t i = 0; i < std::min(v.size(), dim); ++i)
        amps[i] = Complex(v[i] / norm, 0);
    for (size_t i = 0; i < dim && i < reg.size(); ++i)
        reg[i] = amps[i];
}

double QMeans::swapTestSimilarity(const std::vector<double>& a,
                                   const std::vector<double>& b,
                                   int shots) {
    double dot = 0, na = 0, nb = 0;
    size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) { dot += a[i]*b[i]; na += a[i]*a[i]; nb += b[i]*b[i]; }
    double denom = std::sqrt(na * nb);
    double cosine = (denom < 1e-12) ? 0.0 : dot / denom;
    std::uniform_real_distribution<double> noise(
        -1.0/std::sqrt((double)shots),
         1.0/std::sqrt((double)shots));
    return std::max(-1.0, std::min(1.0, cosine + noise(m_rng)));
}

double QMeans::quantumDistance(const std::vector<double>& a,
                                const std::vector<double>& b) {
    double sim = swapTestSimilarity(a, b);
    return std::sqrt(std::max(0.0, 2.0 * (1.0 - sim)));
}

void QMeans::initCentroidsRandom() {
    std::uniform_int_distribution<size_t> d(0, m_data.size() - 1);
    m_centroids.clear();
    std::vector<size_t> chosen;
    while ((int)m_centroids.size() < m_k) {
        size_t idx = d(m_rng);
        if (std::find(chosen.begin(), chosen.end(), idx) == chosen.end()) {
            m_centroids.push_back(m_data[idx]);
            chosen.push_back(idx);
        }
    }
}

void QMeans::assignLabels() {
    for (size_t i = 0; i < m_data.size(); ++i) {
        double best = std::numeric_limits<double>::max();
        int label = 0;
        for (int c = 0; c < m_k; ++c) {
            double dist = quantumDistance(m_data[i], m_centroids[c]);
            if (dist < best) { best = dist; label = c; }
        }
        m_labels[i] = label;
    }
}

void QMeans::updateCentroids() {
    int dim = (int)m_data[0].size();
    std::vector<std::vector<double>> newC(m_k, std::vector<double>(dim, 0.0));
    std::vector<int> counts(m_k, 0);
    for (size_t i = 0; i < m_data.size(); ++i) {
        int c = m_labels[i];
        for (int d = 0; d < dim; ++d) newC[c][d] += m_data[i][d];
        counts[c]++;
    }
    for (int c = 0; c < m_k; ++c)
        if (counts[c] > 0)
            for (int d = 0; d < dim; ++d) newC[c][d] /= counts[c];
    m_centroids = newC;
}

bool QMeans::converged(const std::vector<std::vector<double>>& oldC) const {
    for (int c = 0; c < m_k; ++c) {
        double dist = 0.0;
        for (size_t d = 0; d < m_centroids[c].size(); ++d) {
            double diff = m_centroids[c][d] - oldC[c][d];
            dist += diff * diff;
        }
        if (std::sqrt(dist) > 1e-6) return false;
    }
    return true;
}

void QMeans::run() {
    if (m_data.empty()) throw std::runtime_error("No data set. Call setData() first.");
    if (m_centroids.empty()) initCentroidsRandom();
    m_itersRun = 0;
    for (int iter = 0; iter < 100; ++iter) {
        auto oldC = m_centroids;
        assignLabels();
        updateCentroids();
        m_itersRun++;
        if (converged(oldC)) break;
    }
}

} // namespace QSE