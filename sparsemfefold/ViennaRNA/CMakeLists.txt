# Mateo Gray 2020

cmake_minimum_required(VERSION 3.8)
project(Vienna)

# enable c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("${CMAKE_SOURCE_DIR}")

set(constraints_SOURCE
    constraints/constraints.c
    constraints/hard.c
    constraints/SHAPE.c
    constraints/soft.c
) 
set(datastructures_SOURCE
    datastructures/basic_datastructures.c
    datastructures/char_stream.c
    datastructures/lists.c
)
set(io_SOURCE
    io/file_formats.c
    io/io_utils.c
)
set(landscape_SOURCE
    landscape/move.c
)
set(loops_SOURCE
    loops/external_pf.c
    loops/external.c
    loops/hairpin.c
    loops/internal.c
    loops/multibranch.c
)
set(params_SOURCE
    params/default.c
    params/params.c
    params/io.c
)
set(utils_SOURCE
    utils/cpu.c
    utils/higher_order_functions.c
    utils/string_utils.c
    utils/structure_utils.c
    utils/utils.c
)
set(vienna_SOURCE
    ${constraints_SOURCE} ${datastructures_SOURCE}
    ${io_SOURCE} ${landscape_SOURCE}
    ${loops_SOURCE} ${params_SOURCE}
    ${utils_SOURCE}
    alphabet.c
    boltzmann_sampling.c
    centroid.c
    cofold.c
    commands.c
    dp_matrices.c
    equilibrium_probs.c
    eval.c
    fold_compound.c
    fold.c
    gquad.c
    grammar.c
    MEA.c
    mfe.c
    mm.c
    model.h model.c
    part_func.c
    ribo.c
    sequence.c
    subopt.c
    unstructured_domains.c  
)

add_library(RNA ${vienna_SOURCE})



