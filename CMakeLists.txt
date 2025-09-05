cmake_minimum_required(VERSION 3.10)
project(QRDecompositionTests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


enable_testing()


find_package(GTest REQUIRED)


add_library(qr_decomposition STATIC qr_decomposition.cpp)


add_executable(unit_tests qr_decomposition_test.cpp)
target_link_libraries(unit_tests GTest::GTest GTest::Main qr_decomposition)
add_test(NAME UnitTests COMMAND unit_tests)


add_executable(integration_tests integration_test.cpp)
target_link_libraries(integration_tests qr_decomposition)
add_test(NAME IntegrationTests COMMAND integration_tests)