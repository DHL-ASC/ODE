cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)


# use linalg from NGSolve, otherwise comment it out
# find_package(NGSolve CONFIG REQUIRED
#   HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
# )


include_directories(src external_dependencies/DHL-LinAlg/src external_dependencies/DHL-LinAlg/external_dependencies/HPC/src external_dependencies/DHL-LinAlg/external_dependencies/HPC/concurrentqueue external_dependencies/DHL-LinAlg/external_dependencies/HPC/taskflow)
link_directories(src external_dependencies/DHL-LinAlg/src external_dependencies/DHL-LinAlg/external_dependencies/HPC/src external_dependencies/DHL-LinAlg/external_dependencies/HPC/concurrentqueue external_dependencies/DHL-LinAlg/external_dependencies/HPC/taskflow)

add_subdirectory (src)

add_executable (test_ode demos/test_ode.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/parallelcomputingtf.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/timer.cc)
# target_link_libraries (test_ode PUBLIC bla)

add_executable (test_newmark demos/test_newmark.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/parallelcomputingtf.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/timer.cc)
# target_link_libraries (test_newmark PUBLIC bla)

add_executable (test_alpha demos/test_alpha.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/parallelcomputingtf.cc external_dependencies/DHL-LinAlg/external_dependencies/HPC/src/timer.cc)
# target_link_libraries (test_alpha PUBLIC bla)

add_subdirectory (mass_spring)
