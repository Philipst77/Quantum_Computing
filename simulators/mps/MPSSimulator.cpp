#include "MPSSimulator.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <numeric>

namespace QSE {

static const double INV_SQRT2 = 1.0 / std::sqrt(2.0);

MPSSimulator::Tensor3 MPSSimulator::makeTensor(int d, int chiL, int chiR) const {
    return Tensor3(d, std::vector<std::vector<Complex>>(chiL,
                   std::vector<Complex>(chiR, Complex(0, 0))));
}

MPSSimulator::MPSSimulator(int numQubits, int maxBondDim)
    : m_n(numQubits), m_chi(maxBondDim)
{
    reset();
}

void MPSSimulator::reset() {
    m_tensors.clear();
    m_sv.clear();
    for (int i = 0; i < m_n; ++i) {
        auto T = makeTensor(2, 1, 1);
        T[0][0][0] = Complex(1, 0);
        T[1][0][0] = Complex(0, 0);
        m_tensors.push_back(T);
        m_sv.push_back({1.0});
    }
}

std::vector<int> MPSSimulator::bondDimensions() const {
    std::vector<int> dims;
    for (auto& s : m_sv) dims.push_back((int)s.size());
    return dims;
}

void MPSSimulator::applySingleSiteGate(int q,
    const std::vector<std::vector<Complex>>& gate)
{
    auto& T = m_tensors[q];
    int chiL = T[0].size(), chiR = T[0][0].size();
    auto newT = makeTensor(2, chiL, chiR);
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int l = 0; l < chiL; ++l)
                for (int r = 0; r < chiR; ++r)
                    newT[a][l][r] += gate[a][b] * T[b][l][r];
    m_tensors[q] = newT;
}

void MPSSimulator::applyH(int q) {
    applySingleSiteGate(q, {{INV_SQRT2, INV_SQRT2},{INV_SQRT2,-INV_SQRT2}});
}
void MPSSimulator::applyX(int q) {
    applySingleSiteGate(q, {{0,1},{1,0}});
}
void MPSSimulator::applyY(int q) {
    applySingleSiteGate(q, {{0,Complex(0,-1)},{Complex(0,1),0}});
}
void MPSSimulator::applyZ(int q) {
    applySingleSiteGate(q, {{1,0},{0,-1}});
}
void MPSSimulator::applyS(int q) {
    applySingleSiteGate(q, {{1,0},{0,Complex(0,1)}});
}
void MPSSimulator::applyT(int q) {
    applySingleSiteGate(q, {{1,0},{0,std::exp(Complex(0,M_PI/4.0))}});
}
void MPSSimulator::applyRz(int q, double theta) {
    applySingleSiteGate(q, {{std::exp(Complex(0,-theta/2)),0},
                             {0,std::exp(Complex(0,theta/2))}});
}
void MPSSimulator::applyRy(int q, double theta) {
    double c = std::cos(theta/2), s = std::sin(theta/2);
    applySingleSiteGate(q, {{Complex(c,0),Complex(-s,0)},
                             {Complex(s,0),Complex(c,0)}});
}

MPSSimulator::SVDResult MPSSimulator::thinSVD(
    const std::vector<std::vector<Complex>>& M, int maxRank) const
{
    int rows = M.size(), cols = M[0].size();
    int rank = std::min({rows, cols, maxRank});
    SVDResult res;
    res.U.assign(rows, std::vector<Complex>(rank, 0));
    res.S.assign(rank, 0);
    res.Vt.assign(rank, std::vector<Complex>(cols, 0));

    std::vector<std::vector<Complex>> MdagM(cols, std::vector<Complex>(cols, 0));
    for (int i = 0; i < cols; ++i)
        for (int j = 0; j < cols; ++j)
            for (int k = 0; k < rows; ++k)
                MdagM[i][j] += std::conj(M[k][i]) * M[k][j];

    std::vector<Complex> v(cols, 1.0/std::sqrt((double)cols));
    for (int sv = 0; sv < rank; ++sv) {
        for (int iter = 0; iter < 20; ++iter) {
            std::vector<Complex> Mv(cols, 0);
            for (int i = 0; i < cols; ++i)
                for (int j = 0; j < cols; ++j)
                    Mv[i] += MdagM[i][j] * v[j];
            double norm = 0;
            for (auto& x : Mv) norm += std::norm(x);
            norm = std::sqrt(norm);
            if (norm < 1e-14) break;
            for (auto& x : Mv) x /= norm;
            v = Mv;
        }
        std::vector<Complex> Mv(rows, 0);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                Mv[i] += M[i][j] * v[j];
        double sigma = 0;
        for (auto& x : Mv) sigma += std::norm(x);
        sigma = std::sqrt(sigma);
        res.S[sv] = sigma;
        if (sigma < 1e-14) break;
        for (int i = 0; i < rows; ++i) res.U[i][sv] = Mv[i] / sigma;
        for (int j = 0; j < cols; ++j) res.Vt[sv][j] = std::conj(v[j]);
        for (int i = 0; i < cols; ++i)
            for (int j = 0; j < cols; ++j)
                MdagM[i][j] -= res.S[sv] * res.S[sv] * std::conj(v[i]) * v[j];
        if (sv + 1 < rank) {
            v.assign(cols, 0);
            v[(sv+1) % cols] = 1.0;
            for (int k = 0; k <= sv; ++k) {
                Complex dot = 0;
                for (int j = 0; j < cols; ++j) dot += std::conj(res.Vt[k][j]) * v[j];
                for (int j = 0; j < cols; ++j) v[j] -= dot * res.Vt[k][j];
            }
            double vn = 0;
            for (auto& x : v) vn += std::norm(x);
            vn = std::sqrt(vn);
            if (vn > 1e-14) for (auto& x : v) x /= vn;
        }
    }
    return res;
}

void MPSSimulator::svdAndTruncate(int bond) {
    auto& TL = m_tensors[bond];
    auto& TR = m_tensors[bond + 1];
    int dL = 2, dR = 2;
    int chiL = TL[0].size(), chiM_old = TL[0][0].size(), chiR = TR[0][0].size();

    int rowDim = dL * chiL, colDim = dR * chiR;
    std::vector<std::vector<Complex>> Theta(rowDim, std::vector<Complex>(colDim, 0));
    for (int sL = 0; sL < dL; ++sL)
        for (int sR = 0; sR < dR; ++sR)
            for (int l = 0; l < chiL; ++l)
                for (int m = 0; m < chiM_old; ++m)
                    for (int r = 0; r < chiR; ++r)
                        Theta[sL*chiL+l][sR*chiR+r] += TL[sL][l][m] * TR[sR][m][r];

    auto svd = thinSVD(Theta, m_chi);
    int newBond = 0;
    for (int i = 0; i < (int)svd.S.size(); ++i)
        if (svd.S[i] > 1e-12) ++newBond;
    newBond = std::max(1, std::min(newBond, m_chi));

    m_sv[bond].resize(newBond);
    for (int i = 0; i < newBond; ++i) m_sv[bond][i] = svd.S[i];

    auto newTL = makeTensor(dL, chiL, newBond);
    auto newTR = makeTensor(dR, newBond, chiR);
    for (int sL = 0; sL < dL; ++sL)
        for (int l = 0; l < chiL; ++l)
            for (int m = 0; m < newBond; ++m)
                newTL[sL][l][m] = svd.U[sL*chiL+l][m];
    for (int sR = 0; sR < dR; ++sR)
        for (int m = 0; m < newBond; ++m)
            for (int r = 0; r < chiR; ++r)
                newTR[sR][m][r] = svd.S[m] * svd.Vt[m][sR*chiR+r];

    m_tensors[bond]     = newTL;
    m_tensors[bond + 1] = newTR;
}

void MPSSimulator::applyTwoSiteGate(int q,
    const std::vector<std::vector<Complex>>& gate)
{
    auto& TL = m_tensors[q];
    auto& TR = m_tensors[q + 1];
    int chiL = TL[0].size(), chiM = TL[0][0].size(), chiR = TR[0][0].size();

    auto newTL = makeTensor(2, chiL, chiM);
    auto newTR = makeTensor(2, chiM, chiR);
    auto Theta = makeTensor(2, chiL * 2, chiR);

    for (int sLp = 0; sLp < 2; ++sLp)
        for (int sRp = 0; sRp < 2; ++sRp)
            for (int sL = 0; sL < 2; ++sL)
                for (int sR = 0; sR < 2; ++sR) {
                    Complex g = gate[sLp*2+sRp][sL*2+sR];
                    if (std::abs(g) < 1e-14) continue;
                    for (int l = 0; l < chiL; ++l)
                        for (int m = 0; m < chiM; ++m)
                            for (int r = 0; r < chiR; ++r)
                                Theta[sLp][l*2+sRp][r] += g * TL[sL][l][m] * TR[sR][m][r];
                }

    for (int sLp = 0; sLp < 2; ++sLp)
        for (int l = 0; l < chiL; ++l)
            for (int sRp = 0; sRp < 2; ++sRp)
                for (int r = 0; r < chiR; ++r)
                    newTL[sLp][l][0] += Theta[sLp][l*2+sRp][r];

    m_tensors[q]     = newTL;
    m_tensors[q + 1] = newTR;
    svdAndTruncate(q);
}

void MPSSimulator::applyCNOT(int ctrl, int tgt) {
    if (std::abs(ctrl - tgt) != 1)
        throw std::invalid_argument("applyCNOT requires adjacent qubits.");
    int q = std::min(ctrl, tgt);
    std::vector<std::vector<Complex>> cnot(4, std::vector<Complex>(4, 0));
    if (ctrl < tgt)
        cnot[0][0] = cnot[1][1] = cnot[2][3] = cnot[3][2] = 1;
    else
        cnot[0][0] = cnot[1][3] = cnot[2][2] = cnot[3][1] = 1;
    applyTwoSiteGate(q, cnot);
}

void MPSSimulator::applyTwoSiteCNOT(int ctrl, int tgt) {
    int lo = std::min(ctrl, tgt), hi = std::max(ctrl, tgt);
    for (int i = lo; i < hi - 1; ++i) applyCNOT(i, i+1);
    applyCNOT(hi-1, hi);
    for (int i = hi-2; i >= lo; --i) applyCNOT(i, i+1);
}

double MPSSimulator::probability(int qubit, int outcome) const {
    const auto& T = m_tensors[qubit];
    int chiL = T[0].size(), chiR = T[0][0].size();
    double p = 0;
    for (int l = 0; l < chiL; ++l)
        for (int r = 0; r < chiR; ++r)
            p += std::norm(T[outcome][l][r]);
    if (qubit > 0)
        for (auto& s : m_sv[qubit-1]) p *= s * s;
    return std::min(1.0, std::max(0.0, p));
}

int MPSSimulator::measure(int qubit, std::mt19937& rng) {
    double p1 = probability(qubit, 1);
    std::uniform_real_distribution<double> d(0, 1);
    int outcome = (d(rng) < p1) ? 1 : 0;

    auto& T = m_tensors[qubit];
    int chiL = T[0].size(), chiR = T[0][0].size();
    for (int l = 0; l < chiL; ++l)
        for (int r = 0; r < chiR; ++r)
            T[1-outcome][l][r] = 0;

    double norm = 0;
    for (int l = 0; l < chiL; ++l)
        for (int r = 0; r < chiR; ++r)
            norm += std::norm(T[outcome][l][r]);
    norm = std::sqrt(norm);
    if (norm > 1e-14)
        for (int l = 0; l < chiL; ++l)
            for (int r = 0; r < chiR; ++r)
                T[outcome][l][r] /= norm;
    return outcome;
}

std::vector<int> MPSSimulator::measureAll(std::mt19937& rng) {
    std::vector<int> results(m_n);
    for (int q = 0; q < m_n; ++q) results[q] = measure(q, rng);
    return results;
}

std::vector<Complex> MPSSimulator::toStateVector() const {
    if (m_n > 20) throw std::runtime_error("toStateVector: too many qubits");
    size_t dim = 1ULL << m_n;
    std::vector<Complex> sv(dim, 0);
    for (size_t i = 0; i < dim; ++i) {
        std::vector<Complex> vec = {1.0};
        for (int q = 0; q < m_n; ++q) {
            int bit = (i >> q) & 1;
            const auto& T = m_tensors[q];
            int chiL = T[0].size(), chiR = T[0][0].size();
            std::vector<Complex> newVec(chiR, 0);
            for (int l = 0; l < chiL && l < (int)vec.size(); ++l)
                for (int r = 0; r < chiR; ++r)
                    newVec[r] += vec[l] * T[bit][l][r];
            vec = newVec;
        }
        for (auto& v : vec) sv[i] += v;
    }
    return sv;
}

} // namespace QSE