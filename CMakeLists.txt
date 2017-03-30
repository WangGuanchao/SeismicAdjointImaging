cmake_minimum_required(VERSION 3.0)
project(SeismicAdjointImaging)
file(MAKE_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/bin)

#! Output path
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/bin)

set(CMAKE_CXX_COMPILER "/opt/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64/icpc")
set(CMAKE_C_COMPILER "/opt/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64/icc")

#! Library
set(SJINC lib/sjinc.h)
set(SJLIB lib/sjabc.c lib/sjfile.c lib/sjmalloc.c lib/sjmath.c lib/sjwave.c)
add_library(sj ${SJLIB} ${SGLIB})

#! RUNTIME
find_library(sj ${SeismicAdjointImaging_BINARY_DIR}/lib)

#! test
set(SOURCEFILE test/testsj.c)
add_executable(testsj ${SOURCEFILE} ${SJINC})
target_link_libraries(testsj sj m)

#! sjbin2su
set(SOURCEFILE src/sjbin2su.c)
add_executable(sjbin2su ${SOURCEFILE} ${SJINC})
target_link_libraries(sjbin2su sj m)
#! sjdisplaysurvey2d
set(SOURCEFILE src/sjdisplaysurvey2d.c)
add_executable(sjdisplaysurvey2d ${SOURCEFILE} ${SJINC})
target_link_libraries(sjdisplaysurvey2d sj m)
#! sjsurvey2d
set(SOURCEFILE src/sjsurvey2d.c)
add_executable(sjsurvey2d ${SOURCEFILE} ${SJINC})
target_link_libraries(sjsurvey2d sj m)

#--------------------------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------------------------

set(CMAKE_CXX_COMPILER "/opt/intel/compilers_and_libraries_2017.1.132/linux/mpi/intel64/bin/mpiicpc")
set(CMAKE_C_COMPILER "/opt/intel/compilers_and_libraries_2017.1.132/linux/mpi/intel64/bin/mpiicc")

find_package ( MPI )
if ( MPI_FOUND )
    #! Library
    set(SJINC lib/sjinc.h)
    set(SJLIB lib/sjabc.c lib/sjfile.c lib/sjmalloc.c lib/sjmath.c lib/sjwave.c)
    add_library(mpisj ${SJINC} ${SJLIB})

    #! RUNTIME
    find_library(mpisj ${SeismicAdjointImaging_BINARY_DIR}/lib)

    #! sjmpiawfd2d
    set(SOURCEFILE mpisrc/sjmpiawfd2d.c)
    add_executable(sjmpiawfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiawfd2d mpisj m)

    #! sjmpiarti2d
    set(SOURCEFILE mpisrc/sjmpiarti2d.c mpisrc/sjimage.h)
    add_executable(sjmpiarti2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiarti2d mpisj m)

    #! sjmpiafwi2d
    set(SOURCEFILE mpisrc/sjmpiafwi2d.c mpisrc/sjimage.h)
    add_executable(sjmpiafwi2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiafwi2d mpisj m)

    #! sjmpiartm2d
    set(SOURCEFILE mpisrc/sjmpiartm2d.c mpisrc/sjimage.h)
    add_executable(sjmpiartm2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiartm2d mpisj m)

    #! sjmpilsartm2d
    set(SOURCEFILE mpisrc/sjmpialsrtm2d.c mpisrc/sjimage.h)
    add_executable(sjmpilsartm2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpilsartm2d mpisj m)

    #! sjmpiewfd2d
    set(SOURCEFILE mpisrc/sjmpiewfd2d.c)
    add_executable(sjmpiewfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiewfd2d mpisj m)

    #! sjmpiewsgfd2d
    set(SOURCEFILE mpisrc/sjmpiewsgfd2d.c)
    add_executable(sjmpiewsgfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiewsgfd2d mpisj m)

    #! sjmpiewvssgfd2d
    set(SOURCEFILE mpisrc/sjmpiewvssgfd2d.c)
    add_executable(sjmpiewvssgfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiewvssgfd2d mpisj m)

endif( MPI_FOUND )