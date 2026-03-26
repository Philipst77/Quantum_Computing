#include "Utils.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>
#include <stdexcept>

namespace QSE {
namespace Utils {

// ── Math helpers ─────────────────────────────────────────────

Complex innerProduct(const QubitRegister& a, const QubitRegister& b) {
    if (a.size() != b.size())
        throw std::invalid_argument("Register sizes must match for inner product");
    Complex result{0, 0};
    for (size_t i = 0; i < a.size(); ++i)
        result += std::conj(a[i]) * b[i];
    return result;
}

double fidelity(const QubitRegister& a, const QubitRegister& b) {
    auto ip = innerProduct(a, b);
    return std::norm(ip);
}

bool areApproxEqual(const QubitRegister& a, const QubitRegister& b, double tol) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::abs(a[i] - b[i]) > tol) return false;
    return true;
}

// ── Bitstring utilities ──────────────────────────────────────

std::string toBitstring(size_t value, int numBits) {
    std::string s(numBits, '0');
    for (int i = numBits - 1; i >= 0; --i) {
        s[numBits - 1 - i] = ((value >> i) & 1) ? '1' : '0';
    }
    return s;
}

size_t fromBitstring(const std::string& bs) {
    size_t val = 0;
    for (char c : bs) {
        val = (val << 1) | (c == '1' ? 1 : 0);
    }
    return val;
}

std::vector<int> intToBits(size_t value, int numBits) {
    std::vector<int> bits(numBits);
    for (int i = 0; i < numBits; ++i)
        bits[i] = (value >> i) & 1;
    return bits;
}

size_t bitsToInt(const std::vector<int>& bits) {
    size_t val = 0;
    for (int i = (int)bits.size() - 1; i >= 0; --i)
        val = (val << 1) | (bits[i] & 1);
    return val;
}

int popcount(size_t value) {
    int count = 0;
    while (value) { count += value & 1; value >>= 1; }
    return count;
}

// ── Modular arithmetic ───────────────────────────────────────

long long modpow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

long long gcd(long long a, long long b) {
    while (b) { a %= b; std::swap(a, b); }
    return a;
}

long long modInverse(long long a, long long mod) {
    long long g = gcd(a, mod);
    if (g != 1) throw std::runtime_error("Modular inverse does not exist");
    // Extended Euclidean
    long long x = 0, y = 1, m = mod;
    while (a > 1) {
        long long q = a / m;
        std::swap(m, a %= m);
        std::swap(x, y -= q * x);
    }
    return (y < 0) ? y + mod : y;
}

std::pair<long long, long long> continuedFraction(double x, int maxIter) {
    long long h0 = 1, h1 = 0, k0 = 0, k1 = 1;
    double val = x;
    for (int i = 0; i < maxIter; ++i) {
        long long a = (long long)val;
        long long h2 = a * h0 + h1;
        long long k2 = a * k0 + k1;
        if (std::abs((double)h2 / k2 - x) < 1e-9)
            return {h2, k2};
        h1 = h0; h0 = h2;
        k1 = k0; k0 = k2;
        double remainder = val - a;
        if (std::abs(remainder) < 1e-10) break;
        val = 1.0 / remainder;
    }
    return {h0, k0};
}

// ── Sampling ─────────────────────────────────────────────────

std::map<std::string, int> sampleCounts(QubitRegister& qreg, int shots) {
    std::map<std::string, int> counts;
    auto probs = qreg.probabilities();
    int n = qreg.numQubits();

    static std::mt19937 rng(std::random_device{}());
    std::discrete_distribution<size_t> dist(probs.begin(), probs.end());

    for (int i = 0; i < shots; ++i) {
        size_t outcome = dist(rng);
        counts[toBitstring(outcome, n)]++;
    }
    return counts;
}

std::map<std::string, double> sampleProbabilities(const QubitRegister& qreg) {
    std::map<std::string, double> probs;
    int n = qreg.numQubits();
    for (size_t i = 0; i < qreg.size(); ++i) {
        double p = std::norm(qreg[i]);
        if (p > 1e-10)
            probs[toBitstring(i, n)] = p;
    }
    return probs;
}

std::string mostLikelyOutcome(const QubitRegister& qreg) {
    size_t best = 0;
    double bestProb = 0;
    for (size_t i = 0; i < qreg.size(); ++i) {
        double p = std::norm(qreg[i]);
        if (p > bestProb) { bestProb = p; best = i; }
    }
    return toBitstring(best, qreg.numQubits());
}

// ── Printing ─────────────────────────────────────────────────

void printStateVector(const QubitRegister& qreg, double threshold) {
    std::cout << "State vector (" << qreg.numQubits() << " qubits):\n";
    bool any = false;
    for (size_t i = 0; i < qreg.size(); ++i) {
        double mag = std::abs(qreg[i]);
        if (mag < threshold) continue;
        any = true;
        std::cout << "  |" << toBitstring(i, qreg.numQubits()) << ">  "
                  << std::fixed << std::setprecision(5)
                  << qreg[i].real() << (qreg[i].imag() >= 0 ? "+" : "")
                  << qreg[i].imag() << "i  (p=" << std::norm(qreg[i]) << ")\n";
    }
    if (!any) std::cout << "  (all amplitudes below threshold)\n";
}

void printProbabilities(const QubitRegister& qreg, double threshold) {
    std::cout << "Measurement probabilities:\n";
    for (size_t i = 0; i < qreg.size(); ++i) {
        double p = std::norm(qreg[i]);
        if (p < threshold) continue;
        std::cout << "  |" << toBitstring(i, qreg.numQubits()) << ">  "
                  << std::fixed << std::setprecision(6) << p << "\n";
    }
}

void printHistogram(const std::map<std::string, int>& counts) {
    int total = 0;
    for (auto& [k, v] : counts) total += v;
    int maxCount = 0;
    for (auto& [k, v] : counts) maxCount = std::max(maxCount, v);
    int barWidth = 40;

    std::cout << "Measurement histogram:\n";
    for (auto& [bs, cnt] : counts) {
        int bars = (maxCount > 0) ? (cnt * barWidth / maxCount) : 0;
        std::cout << "  |" << bs << ">  "
                  << std::string(bars, '#')
                  << std::string(barWidth - bars, ' ')
                  << "  " << cnt << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * cnt / total) << "%)\n";
    }
}

// ── Formatting ───────────────────────────────────────────────

std::string complexToString(const Complex& c, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision)
        << c.real() << (c.imag() >= 0 ? "+" : "") << c.imag() << "i";
    return oss.str();
}

std::string phaseToString(double phase) {
    constexpr double PI = M_PI;
    auto approx = [&](double a, double b) { return std::abs(a - b) < 1e-9; };
    if (approx(phase, 0))       return "0";
    if (approx(phase, PI))      return "π";
    if (approx(phase, PI/2))    return "π/2";
    if (approx(phase, PI/4))    return "π/4";
    if (approx(phase, PI/8))    return "π/8";
    if (approx(phase, -PI))     return "-π";
    if (approx(phase, -PI/2))   return "-π/2";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << phase;
    return oss.str();
}

} // namespace Utils
} // namespace QSE