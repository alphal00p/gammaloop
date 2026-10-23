cmake_minimum_required(VERSION 3.16)
project(ec_spqr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(OGDF_LIBRARY_TARGETS_ONLY ON CACHE BOOL "" FORCE)
set(OGDF_DISABLE_HELP_MESSAGE ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(OGDF_ARCH "native" CACHE STRING "Target CPU architecture" FORCE)
if(EMSCRIPTEN)
    set(OGDF_MEMORY_MANAGER "POOL_NTS" CACHE STRING "Single-threaded WASM allocator" FORCE)
    add_compile_options(-fwasm-exceptions)
endif()
add_subdirectory("${OGDF_SOURCE_DIR}" ogdf EXCLUDE_FROM_ALL)

if(EMSCRIPTEN)
    find_program(WASM_OPT wasm-opt REQUIRED)
    add_executable(ec-layout plugin.cpp)
    target_link_options(ec-layout PRIVATE
        -fwasm-exceptions -sSTANDALONE_WASM=1 -sALLOW_MEMORY_GROWTH=1
        -sSTACK_SIZE=1048576 --no-entry
        -Wl,--export=decompose -Wl,--export=initialize -Wl,--export-memory)
    set_target_properties(ec-layout PROPERTIES SUFFIX ".raw.wasm")
    # Expected nonplanarity/insertion failures use explicit status. Typst's wasmi
    # has no exception tags; unexpected errors become traps, like Rust panic=abort.
    add_custom_command(TARGET ec-layout POST_BUILD
        COMMAND "${WASM_OPT}" "$<TARGET_FILE:ec-layout>"
            --enable-exception-handling --enable-bulk-memory --enable-bulk-memory-opt
            --enable-mutable-globals --enable-sign-ext --enable-nontrapping-float-to-int
            --strip-eh --vacuum --dce --remove-unused-module-elements -O2
            -o "${CMAKE_CURRENT_BINARY_DIR}/ec-layout.wasm"
        BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/ec-layout.wasm")
    set(EC_TARGET ec-layout)
else()
    add_executable(ec-spqr main.cpp)
    set(EC_TARGET ec-spqr)
endif()
target_include_directories(${EC_TARGET} SYSTEM PRIVATE
    "$<TARGET_PROPERTY:OGDF,INTERFACE_INCLUDE_DIRECTORIES>" "${JSON_INCLUDE_DIR}")
target_link_libraries(${EC_TARGET} PRIVATE OGDF)
target_compile_options(${EC_TARGET} PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(NOT EMSCRIPTEN)
    enable_testing()
    add_executable(ec-layout-test test_layout.cpp)
    target_include_directories(ec-layout-test SYSTEM PRIVATE
        "$<TARGET_PROPERTY:OGDF,INTERFACE_INCLUDE_DIRECTORIES>" "${JSON_INCLUDE_DIR}")
    target_link_libraries(ec-layout-test PRIVATE OGDF)
    target_compile_options(ec-layout-test PRIVATE -Wall -Wextra -Wpedantic -Werror)
    add_test(NAME ec-layout-reference COMMAND ec-layout-test
        "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/constraints.json"
        "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/seeds.json")
endif()
