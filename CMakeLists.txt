cmake_minimum_required (VERSION 3.22)

set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

set (WHEEL_TYPE "2x3x5" CACHE STRING "Wheel.")
set_property(CACHE WHEEL_TYPE PROPERTY STRINGS "2x3x5" "2x3" "2")


project (Eratosthenes VERSION 1.0.0)

# Support for C++20 is required globally.
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

option (ENABLE_OPENMP "It enables OpenMP, needs a compiler support." ON)
if (ENABLE_OPENMP)
	find_package (OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif ()
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Setup the target optimization with the exception of ARMv8.
    if (NOT CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
        add_compile_options(-march=native)
    else ()
        # Use cortex-a53 tuning since it has non-zero LSL instruction cost
        # and it is resulting in the use of barrel shifts in newer GCCs.
        add_compile_options(-mcpu=native -mtune=cortex-a53)
    endif ()

    # More warnings desired.
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Add possibility to use PGO.
    set(PGO_STAGE "None" CACHE STRING "Profile guided optimization stage.")
    set_property(CACHE PGO_STAGE PROPERTY STRINGS "None" "Collect" "Optimize")
    if (PGO_STAGE STREQUAL "Collect")
	    add_compile_options(-fprofile-dir=./pgo -fprofile-generate=./pgo -flto=auto)
    elseif (PGO_STAGE STREQUAL "Optimize")
	    add_compile_options(-fprofile-dir=./pgo -fprofile-use=./pgo -fprofile-correction -flto=auto)
    endif()
endif()

# Fill metainformation about the CPU and the operating system.
execute_process(COMMAND bash -c "LANG=C lscpu --json | jq '.lscpu[]|select(.field == \"Model name:\")' | jq -r '.data' | tr -d '\\n'" OUTPUT_VARIABLE CMAKE_PROCESSOR_NAME)
execute_process(COMMAND bash -c "LANG=C lscpu --json | jq '.lscpu[]|select(.field == \"Thread(s) per core:\")' | jq -r '.data' | tr -d '\\n'" OUTPUT_VARIABLE CMAKE_THREADS_PER_CORE)
cmake_host_system_information(RESULT CMAKE_IS_64BIT QUERY IS_64BIT)
cmake_host_system_information(RESULT CMAKE_TOTAL_PHYSICAL_MEMORY QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT CMAKE_OS_NAME QUERY OS_NAME)
cmake_host_system_information(RESULT CMAKE_OS_RELEASE QUERY OS_RELEASE)
cmake_host_system_information(RESULT CMAKE_OS_VERSION QUERY OS_VERSION)

# Setup the USE variable based on the selected wheel type.
set(USE_WHEEL_2_3_5 0)
set(USE_WHEEL_2_3 0)
set(USE_WHEEL_2 0)

if ("${WHEEL_TYPE}" STREQUAL "2x3x5")
    set(USE_WHEEL_2_3_5 1)
elseif ("${WHEEL_TYPE}" STREQUAL "2x3")
    set(USE_WHEEL_2_3 1)
else()
    set(USE_WHEEL_2 1)
endif()

configure_file ("cmake_vars.h.in" "${PROJECT_BINARY_DIR}/cmake_vars.h")
include_directories ("${PROJECT_BINARY_DIR}")

add_executable (prime_sieve_bench prime_sieve_bench.cpp eratosthenes_sieve.cpp cpu_info.cpp)
add_executable (show_cpu_caches show_cpu_caches.cpp cpu_info.cpp)
add_executable (factorization_wheel_demo factorization_wheel_demo.cpp)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (PGO_STAGE STREQUAL "Collect")
        target_link_libraries(prime_sieve_bench -lgcov)
        target_link_libraries(show_cpu_caches -lgcov)
        target_link_libraries(factorization_wheel_demo -lgcov)
        add_custom_target(profile ./prime_sieve_bench 10000000000 - COMMAND ./show_cpu_caches COMMAND ./factorization_wheel_demo DEPENDS prime_sieve_bench show_cpu_caches factorization_wheel_demo)
    endif ()
endif()

install (TARGETS prime_sieve_bench DESTINATION bin)
install (TARGETS show_cpu_caches DESTINATION bin)

