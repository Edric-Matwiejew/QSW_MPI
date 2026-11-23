cmake_minimum_required(VERSION 3.20)
project(QSW_MPI LANGUAGES C Fortran)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_compile_options("$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fcheck=all>")
add_library(qsw_core STATIC
    native/iso_precisions.f90
    native/sparse.f90
    native/operators.f90
    native/one_norms.f90
    native/expm.f90
)
set_target_properties(qsw_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(qsw_core PUBLIC MPI::MPI_Fortran)

option(QSW_OPENMP "Enable OpenMP threading" OFF)
if(QSW_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(qsw_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE f2py_include
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)

function(add_f2py_module name source wrapper)
    set(signature "${CMAKE_CURRENT_SOURCE_DIR}/native/${name}.pyf")
    set(kind_map "${CMAKE_CURRENT_SOURCE_DIR}/native/.f2py_f2cmap")
    add_custom_command(
        OUTPUT "${name}module.c" "${wrapper}"
        COMMAND "${Python_EXECUTABLE}" -m numpy.f2py "${signature}" --f2cmap "${kind_map}"
        DEPENDS "${signature}" "${kind_map}"
        VERBATIM
    )
    Python_add_library(${name} MODULE WITH_SOABI
        "${source}" "${name}module.c" "${wrapper}" "${f2py_include}/fortranobject.c"
    )
    target_include_directories(${name} PRIVATE "${f2py_include}")
    target_link_libraries(${name} PRIVATE Python::NumPy ${ARGN})
    install(TARGETS ${name} LIBRARY DESTINATION qsw_mpi)
endfunction()

add_f2py_module(foperators native/foperators.f90 foperators-f2pywrappers.f qsw_core)
add_f2py_module(fMPI native/fMPI.f90 fMPI-f2pywrappers.f qsw_core)

set(QSW_PARALLEL_IO AUTO CACHE STRING "Build the parallel HDF5 writer (AUTO, ON or OFF)")
set_property(CACHE QSW_PARALLEL_IO PROPERTY STRINGS AUTO ON OFF)
if(QSW_PARALLEL_IO)
    set(HDF5_PREFER_PARALLEL TRUE)
    if(APPLE)
        # Keep the IO-extension's HDF5 symbols private as h5py may load another HDF5.
        set(HDF5_USE_STATIC_LIBRARIES TRUE)
    endif()
    find_package(HDF5 QUIET COMPONENTS Fortran)
    if(HDF5_FOUND AND (HDF5_IS_PARALLEL OR HDF5_PROVIDES_PARALLEL))
        add_f2py_module(fHDF5 native/parallel_io.f90 fHDF5-f2pywrappers2.f90 MPI::MPI_Fortran HDF5::HDF5)
        if(APPLE)
            target_link_options(fHDF5 PRIVATE "LINKER:-exported_symbol,_PyInit_fHDF5")
        endif()
    elseif(NOT QSW_PARALLEL_IO STREQUAL "AUTO")
        message(FATAL_ERROR "QSW_PARALLEL_IO requires parallel HDF5 with its Fortran interface. Set HDF5_ROOT to its installation prefix.")
    endif()
endif()
