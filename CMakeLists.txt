cmake_minimum_required(VERSION 3.10)

project(adcircdgswem LANGUAGES Fortran)

enable_language(Fortran)

option(BUILD_DEBUG "Build in debug mode" OFF)
option(USE_MPI "Build parallel version of adcircdgswem" OFF)
option(PYTHON_INTERFACE "Build the Python interface using f2py" OFF)

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

#------------------------------------------#
# Python Interface
if(PYTHON_INTERFACE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    add_library(adcircdgswem_shared_obj OBJECT ${SOURCES})
    add_library(adcircdgswem_shared SHARED $<TARGET_OBJECTS:adcircdgswem_shared_obj>)

    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    find_package(F2PY REQUIRED)

    include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

    set(f2py_module_name "pyadcircdgswem")
    set(pped_src_filename "${f2py_module_name}.F90")
    set(fortran_src_file "${CMAKE_SOURCE_DIR}/srcpy/adcircdgswem_python.F95")
    set(generated_pped_file ${CMAKE_CURRENT_BINARY_DIR}/${pped_src_filename})
    set(generated_module_file "${f2py_module_name}\*.so")

    set(pyadcircdgswem_preproc_FLAGS ${FFLAGS})
    add_custom_target(preprocess_src_file ALL
      COMMAND ${CMAKE_Fortran_COMPILER}
              ${pyadcircdgswem_preproc_FLAGS}
              -I${CMAKE_SOURCE_DIR}/srcpy
              -E ${fortran_src_file}
              -o ${generated_pped_file}
      DEPENDS ${fortran_src_file}
      COMMENT "Preprocessing ${fortran_src_file} for use with f2py."
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/srcpy
      VERBATIM
    )
    set(F2PY_EXECUTABLE ${F2PY_EXECUTABLE} --f90exec=${CMAKE_Fortran_COMPILER})

    add_custom_target(${f2py_module_name} ALL
        COMMAND ${F2PY_EXECUTABLE}
                --verbose ${F2PY_FLAGS}
                -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
                -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
                -ladcircdgswem_shared
                -I${CMAKE_Fortran_MODULE_DIRECTORY}
                -c
                -m ${f2py_module_name}
                ${generated_pped_file}
        DEPENDS ${generated_pped_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        VERBATIM
        )
    add_dependencies(${f2py_module_name} adcircdgswem_shared preprocess_src_file)

    set(PROJECT_PYTHON_DIR ${PROJECT_SOURCE_DIR}/pyADCIRCDGSWEM)
    set(pyadcircdgswempathfile "${PROJECT_PYTHON_DIR}/pyadcircdgswem_path.py")
    file(WRITE ${pyadcircdgswempathfile} "#!/usr/bin/env python3\n")
    string(TIMESTAMP pyadcircdgswemtime)
    file(APPEND ${pyadcircdgswempathfile} "# Auto-generated by CMake on ${pyadcircdgswemtime}\n")
    file(APPEND ${pyadcircdgswempathfile} "# This file will automatically be overwritten on every compilation\n")
    file(APPEND ${pyadcircdgswempathfile} "import os\n")
    file(APPEND ${pyadcircdgswempathfile} "import sys\n")
    file(APPEND ${pyadcircdgswempathfile} "sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))\n")

    add_custom_command(TARGET adcircdgswem_shared POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libadcircdgswem_shared.so
                                         ${PROJECT_PYTHON_DIR}
    )
    add_custom_command(TARGET ${f2py_module_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "${f2py_module_name}\*.so"
                                         ${PROJECT_PYTHON_DIR}
    )

endif(PYTHON_INTERFACE)
