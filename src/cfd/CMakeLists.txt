set(SRC
    cfd27.cpp
    mat/csrmat.cpp
    mat/lodmat.cpp
    mat/densemat.cpp
    mat/sparse_matrix_solver.cpp
    mat/algebraic_bc.cpp
    mat/umfpack_solver.cpp
    geom/simplex.cpp
    geom/searcher.cpp
    geom/jacobi.cpp
    geom/polygon.cpp
    grid/vtk.cpp
    grid/i_grid.cpp
    grid/grid1d.cpp
    grid/regular_grid2d.cpp
    grid/unstructured_grid2d.cpp
    grid/cell_finder.cpp
    debug/printer.cpp
    debug/saver.cpp
    debug/tictoc.cpp
    debug/grid_manip.cpp
    fvm/fvm_extended_collocations.cpp
    fvm/fvm_gradient.cpp
    fvm/fvm_dfdn.cpp
)

add_library(
    ${CFD_TARGET}
    STATIC
    ${SRC}
)

set(UMFPACK_LIBRARIES
    /usr/lib/x86_64-linux-gnu/libumfpack.so
    /usr/lib/x86_64-linux-gnu/libamd.so
    /usr/lib/x86_64-linux-gnu/libcholmod.so
    /usr/lib/x86_64-linux-gnu/libcolamd.so
    /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
    /usr/lib/x86_64-linux-gnu/libblas.so
    m
)

target_link_libraries(${CFD_TARGET} PRIVATE ${UMFPACK_LIBRARIES})
