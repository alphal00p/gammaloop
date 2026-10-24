cmake_minimum_required(VERSION 3.16)
project(ec_spqr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(OGDF_LIBRARY_TARGETS_ONLY ON CACHE BOOL "" FORCE)
set(OGDF_DISABLE_HELP_MESSAGE ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(OGDF_ARCH "native" CACHE STRING "Target CPU architecture" FORCE)
add_subdirectory("${OGDF_SOURCE_DIR}" ogdf EXCLUDE_FROM_ALL)

add_executable(ec-spqr main.cpp)
target_include_directories(ec-spqr SYSTEM PRIVATE
    "$<TARGET_PROPERTY:OGDF,INTERFACE_INCLUDE_DIRECTORIES>" "${JSON_INCLUDE_DIR}")
target_link_libraries(ec-spqr PRIVATE OGDF)
target_compile_options(ec-spqr PRIVATE -Wall -Wextra -Wpedantic -Werror)

enable_testing()
add_executable(ec-layout-test test_layout.cpp)
target_include_directories(ec-layout-test SYSTEM PRIVATE
    "$<TARGET_PROPERTY:OGDF,INTERFACE_INCLUDE_DIRECTORIES>" "${JSON_INCLUDE_DIR}")
target_link_libraries(ec-layout-test PRIVATE OGDF)
target_compile_options(ec-layout-test PRIVATE -Wall -Wextra -Wpedantic -Werror)
add_test(NAME ec-layout-reference COMMAND ec-layout-test
    "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/constraints.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/seeds.json")
