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
    set(FFLAGS2 ${FFLAGS2} -DMPI -DHAVE_MPI_MOD)
    if(BUILD_DEBUG)
        set(FFLAGS2 ${FFLAGS2} -DDEBUG_MPI)
    endif(BUILD_DEBUG)
    #ADCIRC's definitions when MPI is enabled
    #set(FFLAGS2 ${FFLAGS2} -DCMPI -DHAVE_MPI_MOD -DCSCA)
endif(USE_MPI)
add_definitions(${FFLAGS2})
set(FFLAGS ${FFLAGS} ${FFLAGS2})

#------------------------------------------------------------------------------#
# Dependencies
#------------------------------------------#
# ADCIRC
message("ADCIRC_ROOT_DIR=${ADCIRC_ROOT_DIR}")
if(NOT ADCIRC_ROOT_DIR)
    message("The Path, ADCIRC_ROOT_DIR, must be set when invoking CMake.")
endif()
set(ADCIRC_LIB_DIR ${ADCIRC_ROOT_DIR}/work)
set(ADCIRC_INC_DIR ${ADCIRC_ROOT_DIR}/work/odircp)
message("ADCIRC_LIB_DIR=${ADCIRC_LIB_DIR}")
message("ADCIRC_INC_DIR=${ADCIRC_INC_DIR}")
# Link
include_directories(${ADCIRC_INC_DIR})
link_directories(${ADCIRC_LIB_DIR})

#------------------------------------------#
# DGSWEM
message("DGSWEM_ROOT_DIR=${DGSWEM_ROOT_DIR}")
if(NOT DGSWEM_ROOT_DIR)
    message("The Path, DGSWEM_ROOT_DIR, must be set when invoking CMake.")
endif()
set(DGSWEM_LIB_DIR ${DGSWEM_ROOT_DIR}/work)
set(DGSWEM_INC_DIR ${DGSWEM_ROOT_DIR}/work/odir6)
message("DGSWEM_LIB_DIR=${DGSWEM_LIB_DIR}")
message("DGSWEM_INC_DIR=${DGSWEM_INC_DIR}")
# Link
include_directories(${DGSWEM_INC_DIR})
link_directories(${DGSWEM_LIB_DIR})

set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -fvisibility=hidden )

#------------------------------------------------------------------------------#
# Targets
#------------------------------------------#
# Static Library
add_library(adcircdgswem_static_obj OBJECT ${SOURCES})
add_library(adcircdgswem_static STATIC $<TARGET_OBJECTS:adcircdgswem_static_obj>)
target_link_libraries(adcircdgswem_static adc dgswem)

#------------------------------------------#
# Executable
add_executable(${PROJECT_NAME} "src/adcircdgswem_program.F90") 
target_link_libraries(${PROJECT_NAME} adcircdgswem_static)

