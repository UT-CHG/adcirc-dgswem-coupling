cmake_minimum_required(VERSION 3.10)

project(adcircdgswem LANGUAGES Fortran)

enable_language(Fortran)

option(BUILD_DEBUG "Build in debug mode" OFF)
option(USE_MPI "Build parallel version of adcircdgswem" OFF)

# CMake settings
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

# Source Files
file(GLOB SOURCES "src/*.F95")
# Filter out the main program file to compile the rest into a library.
list(FILTER SOURCES EXCLUDE REGEX ".*src/adcircdgswem_program.F95$")

#------------------------------------------------------------------------------#
# Build Options
#------------------------------------------#
# Mandatory definitionas
add_definitions(-DREAL8 -DLINUX)
set(FFLAGS ${FFLAGS} -DREAL8 -DLINUX)

#------------------------------------------#
# Build - Debug
if(BUILD_DEBUG)
    add_definitions(-g -O0 -DDEBUG)
    set(FFLAGS ${FFLAGS} -g -O0 -DDEBUG)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        add_definitions(-fimplicit-none -Wall -fcheck=all -pedantic -fbacktrace)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        add_definitions(-check all -traceback)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
        # using Clang -- To Do
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "MSVC")
        # using Visual Studio C++ -- To Do
    endif()
else(BUILD_DEBUG)
    set(FFLAGS ${FFLAGS} -Wno-unused-variable)
    set(FFLAGS ${FFLAGS} -Wno-unused-function)
    set(FFLAGS ${FFLAGS} -Wno-unused-parameter)
endif(BUILD_DEBUG)

#------------------------------------------#
# Build - MPI
if(USE_MPI)
    find_package(MPI REQUIRED Fortran)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
    set(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
    add_definitions(-DMPI)
    set(FFLAGS ${FFLAGS} -DMPI)
    if(BUILD_DEBUG)
        add_definitions(-DDEBUG_MPI)
        set(FFLAGS ${FFLAGS} -DDEBUG_MPI)
    endif(BUILD_DEBUG)
    #ADCIRC's definitions when MPI is enabled
    #add_definitions(-DCMPI -DHAVE_MPI_MOD -DCSCA)
    #set(FFLAGS ${FFLAGS} -DCMPI -DHAVE_MPI_MOD -DCSCA)
else()
    # Do nothing for now.
endif(USE_MPI)

#------------------------------------------------------------------------------#
# Targets
#------------------------------------------#
# Static Library
add_library(adcircdgswem_static_obj OBJECT ${SOURCES})
add_library(adcircdgswem_static STATIC $<TARGET_OBJECTS:adcircdgswem_static_obj>)

#------------------------------------------#
# Executable
add_executable(${PROJECT_NAME} "src/adcircdgswem_program.F95") 
target_link_libraries(${PROJECT_NAME} adcircdgswem_static)

