#------------------------------------------------------------------------------#
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new
            directory (called a build directory) and run CMake from there.")
endif()

# guard against bad build-type strings
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

#------------------------------------------------------------------------------#
# Build type
if(NOT CMAKE_BUILD_TYPE MATCHES
        "^(|Debug|ExtraDebug|Release|RelWithDebInfo)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE parameter should be left empty or set
            to Debug (-g -DDEBUG), ExtraDebug (-g -DDEBUG -DDEBUG_TEMPORARY),
            Release, RelWithDebInfo (-g).")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    # Unused for now.
    set(FFLAGS ${FFLAGS} -Wno-unused-variable)
    set(FFLAGS ${FFLAGS} -Wno-unused-function)
    set(FFLAGS ${FFLAGS} -Wno-unused-parameter)
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    set(FFLAGS1 ${FFLAGS1} -g)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(FFLAGS1 ${FFLAGS1} -g -O0 -DDEBUG)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        add_definitions(-fimplicit-none -Wall -fcheck=all -pedantic -fbacktrace)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        add_definitions(-check all -traceback)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
        # using Clang -- To Do
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "MSVC")
        # using Visual Studio C++ -- To Do
    endif()
elseif(CMAKE_BUILD_TYPE MATCHES "ExtraDebug")
    set(FFLAGS1 ${FFLAGS1} -g -O0 -DDEBUG -DDEBUG_TEMPORARY)
endif()
add_definitions(${FFLAGS1})
set(FFLAGS ${FFLAGS} ${FFLAGS1})

