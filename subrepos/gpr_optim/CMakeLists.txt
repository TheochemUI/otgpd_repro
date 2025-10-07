project(gprdimer)

cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -DNDEBUG")

# Messages
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message(STATUS "C compiler: " ${CMAKE_C_COMPILER})
message(STATUS "CXX compiler: " ${CMAKE_CXX_COMPILER})

include_directories(./googletest/googletest/include)
link_directories(./googletest/googletest/build/lib)

if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(STATUS "Building ${PROJECT_NAME} as a binary")
  add_executable(
          ${PROJECT_NAME}
          tests/data_types/FieldTest.cpp
          tests/gpr/auxiliary/DistanceTest.cpp
          tests/gpr/auxiliary/ProblemSetUpTest.cpp
          tests/gpr/auxiliary/GradientTest.cpp
          tests/gpr/covariance_functions/SexpatCFTest.cpp
          tests/gpr/dimer/DimerTest.cpp
          tests/gpr/observation_models/LikGaussianTest.cpp
          tests/gpr/ml/GaussianProcessRegressionTest.cpp
          tests/gpr/ml/SCGTest.cpp
          tests/gpr/prior/PriorGaussianTest.cpp
          tests/gpr/prior/PriorLogUnifTest.cpp
          tests/gpr/prior/PriorSqrttTest.cpp
          tests/gpr/prior/PriorTTest.cpp
          tests/gpr/potentials/EAMPotentialTest.cpp
          tests/gpr/AtomicDimerTest.cpp
          tests/managers/io/FileManagerTest.cpp
          backend/DistributionFunctions.cpp
          gpr/auxiliary/Distance.cpp
          gpr/auxiliary/ProblemSetUp.cpp
          gpr/auxiliary/Gradient.cpp
          gpr/covariance_functions/ConstantCF.cpp
          gpr/covariance_functions/SexpatCF.cpp
          gpr/dimer/Dimer.cpp
          gpr/optimization/LBFGS.cpp
          gpr/ml/GaussianProcessRegression.cpp
          gpr/ml/SCG.cpp
          gpr/observation_models/LikGaussian.cpp
          gpr/prior/PriorGaussian.cpp
          gpr/prior/PriorLogUnif.cpp
          gpr/prior/PriorSqrtt.cpp
          gpr/prior/PriorT.cpp
          gpr/prior/PriorUnif.cpp
          gpr/potentials/EAMPotential.cpp
          gpr/AtomicDimer.cpp
          gpr/AtomicDimerInit.cpp
          managers/io/FileManager.cpp
          main.cpp
  )

  # Copy auxiliary  files for tests
  if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Configuring directory 'tests'")
    file(COPY ${CMAKE_SOURCE_DIR}/tests/reference DESTINATION ${CMAKE_BINARY_DIR}/tests)
    file(COPY ${CMAKE_SOURCE_DIR}/mEAMCUH2 DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY ${CMAKE_SOURCE_DIR}/pot.par DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY ${CMAKE_SOURCE_DIR}/input DESTINATION ${CMAKE_BINARY_DIR})
  endif()
else()
  message(STATUS "Building ${PROJECT_NAME} as a shared library")
  add_library(
          ${PROJECT_NAME} SHARED
          backend/DistributionFunctions.cpp
          gpr/auxiliary/Distance.cpp
          gpr/auxiliary/ProblemSetUp.cpp
          gpr/auxiliary/Gradient.cpp
          gpr/covariance_functions/ConstantCF.cpp
          gpr/covariance_functions/SexpatCF.cpp
          gpr/dimer/Dimer.cpp
          gpr/optimization/LBFGS.cpp
          gpr/ml/GaussianProcessRegression.cpp
          gpr/ml/SCG.cpp
          gpr/observation_models/LikGaussian.cpp
          gpr/prior/PriorGaussian.cpp
          gpr/prior/PriorLogUnif.cpp
          gpr/prior/PriorSqrtt.cpp
          gpr/prior/PriorT.cpp
          gpr/prior/PriorUnif.cpp
          gpr/potentials/EAMPotential.cpp
          gpr/AtomicDimer.cpp
          gpr/AtomicDimerInit.cpp
          managers/io/FileManager.cpp
  )
endif()

target_link_libraries(${PROJECT_NAME} gtest)
target_link_libraries(${PROJECT_NAME} gtest_main)
target_link_libraries(${PROJECT_NAME} pthread)

if(FIND_EIGEN)
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
  include_directories(${EIGEN3_INCLUDE_DIR})
  message(STATUS "EIGEN: " ${EIGEN3_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
endif()

# Install
install(TARGETS ${PROJECT_NAME})
