#include "CliffordSimulator.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <random>

namespace QSE {

CliffordSimulator::CliffordSimulator(int numQubits) : n_(numQubits) {
    tab_.assign(2 * n_, std::vector<int>(2 * n_ + 1, 0));
    for (int i = 0; i < n_; ++i) {
        tab_[i][i]           = 1; // destabiliser: X_i
        tab_[i + n_][i + n_] = 1; // stabiliser:   Z_i
    }
}

void CliffordSimulator::reset() {
    tab_.assign(2 * n_, std::vector<int>(2 * n_ + 1, 0));
    for (int i = 0; i < n_; ++i) {
        tab_[i][i]           = 1;
        tab_[i + n_][i + n_] = 1;
    }
}

int CliffordSimulator::g(int x1, int z1, int x2, int z2) const {
    if (!x1 && !z1) return 0;
    if (x1 && !z1)  return z2 * (1 - 2*x2);
    if (!x1 && z1)  return x2 * (2*z2 - 1);
    return z2 - x2;
}

void CliffordSimulator::rowsum(int h, int i) {
    int phase = 2 * tab_[h][2*n_] + 2 * tab_[i][2*n_];
    for (int j = 0; j < n_; ++j)
        phase += g(tab_[i][j], tab_[i][j+n_],
                   tab_[h][j], tab_[h][j+n_]);
    phase = ((phase % 4) + 4) % 4;
    tab_[h][2*n_] = (phase == 2 || phase == 3) ? 1 : 0;
    for (int j = 0; j < n_; ++j) {
        tab_[h][j]       ^= tab_[i][j];
        tab_[h][j + n_]  ^= tab_[i][j + n_];
    }
}

void CliffordSimulator::h(int q) {
    for (int i = 0; i < 2 * n_; ++i) {
        tab_[i][2*n_] ^= tab_[i][q] & tab_[i][q + n_];
        std::swap(tab_[i][q], tab_[i][q + n_]);
    }
}

void CliffordSimulator::s(int q) {
    for (int i = 0; i < 2 * n_; ++i) {
        tab_[i][2*n_]  ^= tab_[i][q] & tab_[i][q + n_];
        tab_[i][q + n_] ^= tab_[i][q];
    }
}

void CliffordSimulator::sdg(int q) { s(q); s(q); s(q); }

void CliffordSimulator::z(int q)  { s(q); s(q); }

void CliffordSimulator::x(int q)  { h(q); z(q); h(q); }

void CliffordSimulator::y(int q)  { x(q); z(q); }

void CliffordSimulator::cx(int ctrl, int tgt) {
    for (int i = 0; i < 2 * n_; ++i) {
        tab_[i][2*n_]    ^= tab_[i][ctrl] & tab_[i][tgt + n_] &
                            (tab_[i][tgt] ^ tab_[i][ctrl + n_] ^ 1);
        tab_[i][tgt]      ^= tab_[i][ctrl];
        tab_[i][ctrl + n_] ^= tab_[i][tgt + n_];
    }
}

void CliffordSimulator::cz(int q1, int q2) { h(q2); cx(q1, q2); h(q2); }

void CliffordSimulator::cy(int q1, int q2) { sdg(q2); cx(q1, q2); s(q2); }

void CliffordSimulator::swap(int q1, int q2) { cx(q1,q2); cx(q2,q1); cx(q1,q2); }

int CliffordSimulator::measure(int qubit) {
    // Find stabiliser row with X on this qubit
    int p = -1;
    for (int i = n_; i < 2 * n_; ++i)
        if (tab_[i][qubit]) { p = i; break; }

    if (p != -1) {
        for (int i = 0; i < 2 * n_; ++i)
            if (i != p && tab_[i][qubit]) rowsum(i, p);

        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> coin(0, 1);
        int outcome = coin(rng);

        std::fill(tab_[p].begin(), tab_[p].end(), 0);
        tab_[p][qubit + n_] = 1;
        tab_[p][2*n_]       = outcome;

        std::fill(tab_[p - n_].begin(), tab_[p - n_].end(), 0);
        tab_[p - n_][qubit] = 1;
        return outcome;
    }

    // Deterministic
    tab_.push_back(std::vector<int>(2*n_+1, 0));
    tab_.push_back(std::vector<int>(2*n_+1, 0));
    for (int i = 0; i < n_; ++i)
        if (tab_[i][qubit]) rowsum(2*n_, i + n_);
    int outcome = tab_[2*n_][2*n_];
    tab_.pop_back();
    tab_.pop_back();
    return outcome;
}

std::string CliffordSimulator::stabilizers() const {
    std::ostringstream oss;
    for (int i = n_; i < 2 * n_; ++i) {
        oss << (tab_[i][2*n_] ? "-" : "+");
        for (int j = 0; j < n_; ++j) {
            int x = tab_[i][j], z = tab_[i][j + n_];
            if (!x && !z) oss << "I";
            else if (x && !z) oss << "X";
            else if (!x && z) oss << "Z";
            else oss << "Y";
        }
        oss << "\n";
    }
    return oss.str();
}

std::string CliffordSimulator::destabilizers() const {
    std::ostringstream oss;
    for (int i = 0; i < n_; ++i) {
        oss << (tab_[i][2*n_] ? "-" : "+");
        for (int j = 0; j < n_; ++j) {
            int x = tab_[i][j], z = tab_[i][j + n_];
            if (!x && !z) oss << "I";
            else if (x && !z) oss << "X";
            else if (!x && z) oss << "Z";
            else oss << "Y";
        }
        oss << "\n";
    }
    return oss.str();
}

std::string CliffordSimulator::fullTableau() const {
    return "Destabilizers:\n" + destabilizers() +
           "Stabilizers:\n"   + stabilizers();
}

} // namespace QSE