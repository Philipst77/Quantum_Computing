#include <iostream>
#include <vector>
#include <algorithm>
#include "../algorithms/grover/Grover.h"

using namespace QSE;

int main() {
    std::cout << "=== Grover's Search Demo ===\n\n";

    int target    = 5;
    int numQubits = 3;

    std::cout << "Searching " << (1 << numQubits) << " states for target |"
              << target << "> ...\n";

    Grover grover(numQubits,
                  [target](size_t state) { return state == (size_t)target; });

    std::cout << "Optimal iterations: " << Grover::optimalIterations(numQubits) << "\n";

    grover.run();

    std::cout << "Found state : " << grover.result() << "\n";
    std::cout << "Correct     : " << (grover.success() ? "YES" : "NO") << "\n";
    std::cout << "\n" << grover.resultSummary() << "\n";

    // ── Example 2: 4-qubit search with multiple targets ──
    std::cout << "\n--- 4-Qubit Search: targets {3, 11} ---\n";
    std::vector<size_t> targets = {3, 11};
    Grover grover2(4,
        [&targets](size_t state) {
            return std::find(targets.begin(), targets.end(), state) != targets.end();
        },
        Grover::optimalIterations(4, 2));

    grover2.run();
    std::cout << "Found state: " << grover2.result() << "\n";

    bool found = std::find(targets.begin(), targets.end(), grover2.result()) != targets.end();
    std::cout << "In target set: " << (found ? "YES" : "NO") << "\n";

    return 0;
}