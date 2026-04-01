#include "QubitRegister.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <numeric>
#include <cassert>

namespace QSE {

QubitRegister::QubitRegister(int numQubits)
    : numQubits_(numQubits),
      state_(1ULL << numQubits, {0.0, 0.0}),
      rng_(std::random_device{}())
{
    if (numQubits <= 0)
        throw std::invalid_argument("Number of qubits must be positive");
    if (numQubits > 28)
        throw std::invalid_argument("Too many qubits (max 28 for statevector simulation)");
    state_[0] = {1.0, 0.0};
}

void QubitRegister::reset() {
    std::fill(state_.begin(), state_.end(), Complex{0.0, 0.0});
    state_[0] = {1.0, 0.0};
}

void QubitRegister::setState(const StateVector& sv) {
    if (sv.size() != stateSize())
        throw std::invalid_argument("State vector size mismatch");
    state_ = sv;
}

void QubitRegister::setComputationalBasis(size_t basisState) {
    if (basisState >= stateSize())
        throw std::out_of_range("Basis state index out of range");
    std::fill(state_.begin(), state_.end(), Complex{0.0, 0.0});
    state_[basisState] = {1.0, 0.0};
}

void QubitRegister::validateQubit(int qubit) const {
    if (qubit < 0 || qubit >= numQubits_)
        throw std::out_of_range("Qubit index out of range");
}

void QubitRegister::applyGate(int qubit,
    const std::vector<std::vector<Complex>>& matrix) {

    validateQubit(qubit);

    size_t n = stateSize();
    size_t bit = 1ULL << qubit;

    for (size_t i = 0; i < n; i += bit * 2) {
        for (size_t j = i; j < i + bit; ++j) {
            Complex a = state_[j];
            Complex b = state_[j + bit];

            state_[j]       = matrix[0][0] * a + matrix[0][1] * b;
            state_[j + bit] = matrix[1][0] * a + matrix[1][1] * b;
        }
    }
}

void QubitRegister::applyControlledGate(int control, int target,
    const std::vector<std::vector<Complex>>& matrix) {

    validateQubit(control);
    validateQubit(target);

    size_t n = stateSize();
    size_t cBit = 1ULL << control;
    size_t tBit = 1ULL << target;

    for (size_t i = 0; i < n; ++i) {
        // Only act when control = 1 and target = 0
        // Process each pair once from the target=0 side
        if ((i & cBit) == 0) continue;
        if ((i & tBit) != 0) continue;

        size_t j = i | tBit;  // partner state with target = 1

        Complex a = state_[i];
        Complex b = state_[j];

        state_[i] = matrix[0][0] * a + matrix[0][1] * b;
        state_[j] = matrix[1][0] * a + matrix[1][1] * b;
    }
}

void QubitRegister::applyToffoliGate(int control1, int control2, int target) {
    validateQubit(control1);
    validateQubit(control2);
    validateQubit(target);

    size_t n = stateSize();
    size_t c1Bit = 1ULL << control1;
    size_t c2Bit = 1ULL << control2;
    size_t tBit  = 1ULL << target;

    for (size_t i = 0; i < n; ++i) {
        // Both controls must be 1, target must be 0
        // Process each pair once from the target=0 side
        if (!(i & c1Bit)) continue;
        if (!(i & c2Bit)) continue;
        if (  i & tBit  ) continue;

        size_t j = i | tBit;  // partner state with target = 1
        std::swap(state_[i], state_[j]);
    }
}

int QubitRegister::measure(int qubit) {
    validateQubit(qubit);

    size_t bit = 1ULL << qubit;

    double prob1 = 0.0;
    for (size_t i = 0; i < stateSize(); ++i)
        if (i & bit)
            prob1 += std::norm(state_[i]);

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int outcome = (dist(rng_) < prob1) ? 1 : 0;

    double normFactor = 1.0 / std::sqrt(outcome ? prob1 : (1.0 - prob1));

    for (size_t i = 0; i < stateSize(); ++i) {
        bool isOne = (i & bit) != 0;

        if (isOne != (bool)outcome)
            state_[i] = {0.0, 0.0};
        else
            state_[i] *= normFactor;
    }

    return outcome;
}

std::vector<int> QubitRegister::measureAll() {
    std::vector<double> probs = probabilities();

    std::discrete_distribution<size_t> dist(probs.begin(), probs.end());
    size_t outcome = dist(rng_);

    std::fill(state_.begin(), state_.end(), Complex{0.0, 0.0});
    state_[outcome] = {1.0, 0.0};

    std::vector<int> bits(numQubits_);
    for (int i = 0; i < numQubits_; ++i)
        bits[i] = (outcome >> i) & 1;

    return bits;
}

double QubitRegister::probability(size_t basisState) const {
    return std::norm(state_[basisState]);
}

std::vector<double> QubitRegister::probabilities() const {
    std::vector<double> probs(stateSize());
    for (size_t i = 0; i < stateSize(); ++i)
        probs[i] = std::norm(state_[i]);
    return probs;
}

void QubitRegister::normalize() {
    double n = norm();
    if (n < 1e-15)
        throw std::runtime_error("Cannot normalize zero-norm state");

    for (auto& amp : state_)
        amp /= n;
}

double QubitRegister::norm() const {
    double n2 = 0.0;
    for (const auto& amp : state_)
        n2 += std::norm(amp);
    return std::sqrt(n2);
}

bool QubitRegister::isNormalized(double tol) const {
    return std::abs(norm() - 1.0) < tol;
}

std::string QubitRegister::toString() const {
    std::ostringstream oss;
    bool first = true;

    for (size_t i = 0; i < stateSize(); ++i) {
        if (std::abs(state_[i]) < 1e-10) continue;

        if (!first) oss << " + ";

        oss << "(" << state_[i].real();
        if (state_[i].imag() >= 0) oss << "+";
        oss << state_[i].imag() << "i)|";

        for (int q = numQubits_ - 1; q >= 0; --q)
            oss << ((i >> q) & 1);

        oss << ">";
        first = false;
    }

    return oss.str();
}

void QubitRegister::print() const {
    std::cout << toString() << std::endl;
}

QubitRegister QubitRegister::tensorProduct(const QubitRegister& a,
                                           const QubitRegister& b) {

    int totalQubits = a.numQubits_ + b.numQubits_;
    QubitRegister result(totalQubits);

    size_t sizeA = a.stateSize();
    size_t sizeB = b.stateSize();

    for (size_t i = 0; i < sizeA; ++i)
        for (size_t j = 0; j < sizeB; ++j)
            result.state_[i * sizeB + j] = a.state_[i] * b.state_[j];

    return result;
}

} // namespace QSE