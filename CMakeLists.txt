cmake_minimum_required(VERSION 3.0.0)
project(Astraea VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()
add_executable(Astraea main.cpp Astraea.cpp Astraea.h TaylorConstants.h lib/Quartic/quartic.cpp lib/Quartic/quartic.h)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
