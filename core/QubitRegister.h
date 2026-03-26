#pragma once
#include <vector>
#include <complex>
#include <string>
#include <random>
#include <stdexcept>
#include <functional>

namespace QSE {

using Complex = std::complex<double>;
using StateVector = std::vector<Complex>;

class QubitRegister {
public:
    explicit QubitRegister(int numQubits);
    QubitRegister(const QubitRegister&) = default;
    QubitRegister& operator=(const QubitRegister&) = default;
    QubitRegister(QubitRegister&&) = default;
    QubitRegister& operator=(QubitRegister&&) = default;

    // State access
    int numQubits() const { return numQubits_; }
    size_t size() const { return state_.size(); }
    const StateVector& state() const { return state_; }
    StateVector& state() { return state_; }
    Complex& operator[](size_t idx) { return state_[idx]; }
    const Complex& operator[](size_t idx) const { return state_[idx]; }

    // State initialization
    void reset();
    void setState(const StateVector& sv);
    void setComputationalBasis(size_t basisState);

    // Single-qubit gate application
    void applyGate(int qubit, const std::vector<std::vector<Complex>>& matrix);

    // Two-qubit gate application
    void applyControlledGate(int control, int target, const std::vector<std::vector<Complex>>& matrix);

    // Three-qubit gate application
    void applyToffoliGate(int control1, int control2, int target);

    // Measurement
    int measure(int qubit);
    std::vector<int> measureAll();
    double probability(size_t basisState) const;
    std::vector<double> probabilities() const;

    // Utilities
    void normalize();
    double norm() const;
    bool isNormalized(double tol = 1e-10) const;
    std::string toString() const;
    void print() const;

    // Tensor product
    static QubitRegister tensorProduct(const QubitRegister& a, const QubitRegister& b);

private:
    int numQubits_;
    StateVector state_;
    std::mt19937 rng_;

    size_t stateSize() const { return 1ULL << numQubits_; }
    void validateQubit(int qubit) const;
};

} // namespace QSE