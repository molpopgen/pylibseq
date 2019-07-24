set (CPPTESTS_SOURCES cpptests.cc
    allele_count_matrix.cc)

pybind11_add_module(cpptests ${CPPTESTS_SOURCES})
set_target_properties(cpptests PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
