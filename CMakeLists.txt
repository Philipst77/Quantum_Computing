cmake_minimum_required(VERSION 3.16)
project(QuantumStateEngine VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ────────────────────────────────────────────────────────────
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-O2 -Wall -Wextra -Wpedantic)
endif()

# ── Core library ──────────────────────────────────────────────────────────────
add_library(qse_core STATIC
    core/QubitRegister.cpp
    core/QuantumGate.cpp
    core/QuantumAlgorithm.cpp
    core/Utils.cpp
)
target_include_directories(qse_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ── Algorithm library ─────────────────────────────────────────────────────────
add_library(qse_algorithms STATIC
    algorithms/grover/Grover.cpp
    algorithms/qft/QFT.cpp
    algorithms/teleportation/Teleportation.cpp
    algorithms/cryptography/BB84.cpp
    algorithms/shor/Shor.cpp
    algorithms/qml/QMeans.cpp
)
target_link_libraries(qse_algorithms PUBLIC qse_core)

# ── Error correction library ──────────────────────────────────────────────────
add_library(qse_error_correction STATIC
    error_correction/ShorCode.cpp
    error_correction/ThreeQubitFlip.cpp
)
target_link_libraries(qse_error_correction PUBLIC qse_core)

# ── Optimization library ──────────────────────────────────────────────────────
add_library(qse_optimization STATIC
    optimization/QAOA.cpp
    optimization/VQE.cpp
)
target_link_libraries(qse_optimization PUBLIC qse_core)

# ── Simulator libraries ───────────────────────────────────────────────────────
add_library(qse_simulators STATIC
    simulators/statevector/StatevectorSimulator.cpp
    simulators/clifford/CliffordSimulator.cpp
    simulators/mps/MPSSimulator.cpp
)
target_link_libraries(qse_simulators PUBLIC qse_core)

# ── All libs combined (convenience) ──────────────────────────────────────────
set(QSE_ALL_LIBS
    qse_algorithms
    qse_error_correction
    qse_optimization
    qse_simulators
)

# ── Examples ──────────────────────────────────────────────────────────────────
add_executable(example_bell_state   examples/bell_state.cpp)
add_executable(example_grover       examples/grover.cpp)
add_executable(example_teleportation examples/teleportation.cpp)

target_link_libraries(example_bell_state    ${QSE_ALL_LIBS})
target_link_libraries(example_grover        ${QSE_ALL_LIBS})
target_link_libraries(example_teleportation ${QSE_ALL_LIBS})

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()

add_executable(test_core       tests/test_core.cpp)
add_executable(test_gates      tests/test_gates.cpp)
add_executable(test_algorithms tests/test_algorithms.cpp)

target_link_libraries(test_core       ${QSE_ALL_LIBS})
target_link_libraries(test_gates      ${QSE_ALL_LIBS})
target_link_libraries(test_algorithms ${QSE_ALL_LIBS})

add_test(NAME CoreTests       COMMAND test_core)
add_test(NAME GateTests       COMMAND test_gates)
add_test(NAME AlgorithmTests  COMMAND test_algorithms)

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS qse_core qse_algorithms qse_error_correction
                qse_optimization qse_simulators
        ARCHIVE DESTINATION lib)

install(DIRECTORY core algorithms error_correction optimization simulators
        DESTINATION include/qse
        FILES_MATCHING PATTERN "*.h")

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "Quantum State Engine (QSE) v${PROJECT_VERSION}")
message(STATUS "  C++ standard : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "  Build with:")
message(STATUS "    mkdir build && cd build")
message(STATUS "    cmake .. && make -j$(nproc)")
message(STATUS "    ctest --output-on-failure")
message(STATUS "")