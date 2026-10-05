set(SRC
    utils/filesystem.cpp
    utils/vecmat.cpp
    cfd27_test.cpp
    infra_test.cpp
    grid_test.cpp
    mat_test.cpp
    geom_test.cpp
    poisson_fdm_test.cpp
    poisson_fvm_test.cpp
    poisson_fvmex_test.cpp
)

add_executable(
    cfd27_test
    ${SRC}
    ${HEADERS}
)
target_link_libraries(cfd27_test PRIVATE ${CFD_TARGET})
target_link_libraries(cfd27_test PRIVATE Catch2::Catch2)
target_include_directories(cfd27_test PRIVATE ${CMAKE_SOURCE_DIR}/src/test)
target_compile_definitions(cfd27_test PUBLIC TEST_DIRECTORY="${CMAKE_SOURCE_DIR}/test_data/")
