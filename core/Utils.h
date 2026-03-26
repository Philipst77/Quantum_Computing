#pragma once
#include "QubitRegister.h"
#include <string>
#include <vector>
#include <map>

namespace QSE {
namespace Utils {

// ── Math helpers ─────────────────────────────────────────────
double fidelity(const QubitRegister& a, const QubitRegister& b);
Complex innerProduct(const QubitRegister& a, const QubitRegister& b);
bool areApproxEqual(const QubitRegister& a, const QubitRegister& b, double tol = 1e-9);

// ── Bitstring utilities ──────────────────────────────────────
std::string toBitstring(size_t value, int numBits);
size_t fromBitstring(const std::string& bs);
std::vector<int> intToBits(size_t value, int numBits);
size_t bitsToInt(const std::vector<int>& bits);
int popcount(size_t value);

// ── Modular arithmetic (for Shor's algorithm) ────────────────
long long modpow(long long base, long long exp, long long mod);
long long modInverse(long long a, long long mod);
long long gcd(long long a, long long b);
std::pair<long long, long long> continuedFraction(double x, int maxIter = 20);

// ── Sampling ─────────────────────────────────────────────────
std::map<std::string, int> sampleCounts(QubitRegister& qreg, int shots);
std::map<std::string, double> sampleProbabilities(const QubitRegister& qreg);
std::string mostLikelyOutcome(const QubitRegister& qreg);

// ── State printing ───────────────────────────────────────────
void printStateVector(const QubitRegister& qreg, double threshold = 1e-10);
void printProbabilities(const QubitRegister& qreg, double threshold = 1e-10);
void printHistogram(const std::map<std::string, int>& counts);

// ── Circuit description helpers ──────────────────────────────
std::string phaseToString(double phase);
std::string complexToString(const Complex& c, int precision = 4);

} // namespace Utils
} // namespace QSE