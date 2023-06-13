project(SparseRNAFolD)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -fconcepts \ -O3")

set(sparseRNAfolD_SOURCE
    src/cmdline.cc
    src/SparseRNAFolD.cc
    src/trace_arrow.cc 
)

set(test_SOURCE
test/Case7.cpp 
test/Test.cpp
)

set(constraints_SOURCE
    src/ViennaRNA/constraints/constraints.c
    src/ViennaRNA/constraints/hard.c
    src/ViennaRNA/constraints/SHAPE.c
    src/ViennaRNA/constraints/soft.c
) 
set(datastructures_SOURCE
    src/ViennaRNA/datastructures/basic_datastructures.c
    src/ViennaRNA/datastructures/char_stream.c
    src/ViennaRNA/datastructures/lists.c
)
set(io_SOURCE
    src/ViennaRNA/io/file_formats.c
    src/ViennaRNA/io/io_utils.c
)
set(landscape_SOURCE
    src/ViennaRNA/landscape/move.c
)
set(loops_SOURCE
    src/ViennaRNA/loops/external_pf.c
    src/ViennaRNA/loops/external.c
    src/ViennaRNA/loops/hairpin.c
    src/ViennaRNA/loops/internal.c
    src/ViennaRNA/loops/multibranch.c
)
set(params_SOURCE
    src/ViennaRNA/params/default.c
    src/ViennaRNA/params/params.c
    src/ViennaRNA/params/io.c
)
set(utils_SOURCE
    src/ViennaRNA/utils/cpu.c
    src/ViennaRNA/utils/higher_order_functions.c
    src/ViennaRNA/utils/string_utils.c
    src/ViennaRNA/utils/structure_utils.c
    src/ViennaRNA/utils/utils.c
)
set(vienna_SOURCE
    ${constraints_SOURCE} ${datastructures_SOURCE}
    ${io_SOURCE} ${landscape_SOURCE}
    ${loops_SOURCE} ${params_SOURCE}
    ${utils_SOURCE}
    src/ViennaRNA/alphabet.c
    src/ViennaRNA/boltzmann_sampling.c
    src/ViennaRNA/centroid.c
    src/ViennaRNA/cofold.c
    src/ViennaRNA/commands.c
    src/ViennaRNA/dp_matrices.c
    src/ViennaRNA/equilibrium_probs.c
    src/ViennaRNA/eval.c
    src/ViennaRNA/fold_compound.c
    src/ViennaRNA/fold.c
    src/ViennaRNA/gquad.c
    src/ViennaRNA/grammar.c
    src/ViennaRNA/MEA.c
    src/ViennaRNA/mfe.c
    src/ViennaRNA/mm.c
    src/ViennaRNA/model.c
    src/ViennaRNA/part_func.c
    src/ViennaRNA/ribo.c
    src/ViennaRNA/sequence.c
    src/ViennaRNA/subopt.c
    src/ViennaRNA/unstructured_domains.c  
)
add_library(RNA STATIC ${vienna_SOURCE})
# configure_file(RNA.pc.in RNA.pc @ONLY)

target_include_directories(RNA PRIVATE .)

add_executable(SparseRNAFolD ${sparseRNAfolD_SOURCE})

add_executable(tests ${test_SOURCE})

target_link_libraries(SparseRNAFolD PRIVATE RNA)

target_link_libraries(tests PRIVATE RNA)


include_directories(src)
