#pragma once
#include "QubitRegister.h"
#include "QuantumGate.h"
#include <iostream>   
#include "Utils.h"
#include <string>
#include <memory>
#include <functional>

namespace QSE {

// Abstract base class for all quantum algorithms
class QuantumAlgorithm {
public:
    explicit QuantumAlgorithm(const std::string& name, int numQubits)
        : name_(name), numQubits_(numQubits), qreg_(numQubits), verbose_(false) {}

    virtual ~QuantumAlgorithm() = default;

    // Core interface — implement in subclasses
    virtual void buildCircuit() = 0;   // Construct the quantum circuit
    virtual void run() = 0;            // Execute and collect results
    virtual std::string resultSummary() const = 0; // Human-readable result

    // Accessors
    const std::string& name() const { return name_; }
    int numQubits() const { return numQubits_; }
    const QubitRegister& qreg() const { return qreg_; }

    void setVerbose(bool v) { verbose_ = v; }
    bool verbose() const { return verbose_; }

    // Reset register to |0...0>
    void reset() { qreg_.reset(); }

    // Convenience: run N shots and return counts
    std::map<std::string, int> sample(int shots = 1024) {
        return Utils::sampleCounts(qreg_, shots);
    }

protected:
    std::string name_;
    int numQubits_;
    QubitRegister qreg_;
    bool verbose_;

    void log(const std::string& msg) const {
        if (verbose_) std::cout << "[" << name_ << "] " << msg << "\n";
    }
};

} // namespace QSE