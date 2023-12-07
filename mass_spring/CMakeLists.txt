add_executable (test_mass_spring mass_spring.cc ../external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/parallelcomputingtf.cc ../external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/timer.cc)
# target_link_libraries (test_mass_spring PUBLIC)


find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(mass_spring bind_mass_spring.cc ../external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/parallelcomputingtf.cc ../external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/timer.cc)
target_link_libraries (mass_spring PUBLIC)
