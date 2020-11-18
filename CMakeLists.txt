cmake_minimum_required(VERSION 3.10)

project(adcircdgswem LANGUAGES Fortran)

enable_language(Fortran)

option(USE_MPI "Build parallel version of adcircdgswem" OFF)

# CMake settings
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

#------------------------------------------------------------------------------#
# Build type
include(BuildType) # For CMAKE_BUILD_TYPE=<options>

if(CMAKE_BUILD_TYPE MATCHES "Debug"
        OR CMAKE_BUILD_TYPE MATCHES "ExtraDebug")
    set(BUILD_DEBUG "ON")
endif()

#------------------------------------------------------------------------------#
# Source Files
file(GLOB SOURCES "src/*.F90")
# Filter out the main program file to compile the rest into a library.
list(FILTER SOURCES EXCLUDE REGEX ".*src/adcircdgswem_program.F90$")

#------------------------------------------------------------------------------#
# Build Options
#------------------------------------------#
# Mandatory definitionas
set(FFLAGS2 ${FFLAGS2} -DREAL8 -DLINUX)

#------------------------------------------#
# Build - MPI
if(USE_MPI)
    find_package(MPI REQUIRED Fortran)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
    set(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
    set(FFLAGS2 ${FFLAGS2} -DMPI)
    if(BUILD_DEBUG)
        set(FFLAGS2 ${FFLAGS2} -DDEBUG_MPI)
    endif(BUILD_DEBUG)
    #ADCIRC's definitions when MPI is enabled
    #set(FFLAGS2 ${FFLAGS2} -DCMPI -DHAVE_MPI_MOD -DCSCA)
endif(USE_MPI)
add_definitions(${FFLAGS2})
set(FFLAGS ${FFLAGS} ${FFLAGS2})

#------------------------------------------------------------------------------#
# Targets
#------------------------------------------#
# Static Library
add_library(adcircdgswem_static_obj OBJECT ${SOURCES})
add_library(adcircdgswem_static STATIC $<TARGET_OBJECTS:adcircdgswem_static_obj>)

#------------------------------------------#
# Executable
add_executable(${PROJECT_NAME} "src/adcircdgswem_program.F90") 
target_link_libraries(${PROJECT_NAME} adcircdgswem_static)

