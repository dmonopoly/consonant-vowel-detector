# CMakeLists for Quixote Cipher.
cmake_minimum_required(VERSION 2.8)

# Options. Turn on with cmake -Dmyvarname=ON.
option(test "Build all tests." OFF)

set(PROJECT_NAME consonant-vowel-detector)
project(${PROJECT_NAME})

# Require c++11 and libc++ for more helpers like initializer_list.
set(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -stdlib=libc++")

# Libraries -> Executables -> Linking.
add_library(basic_helper_lib EMViterbiPackage/BasicHelper.cc)
add_library(notation_lib EMViterbiPackage/Notation.cc)
add_library(notation_helper_lib EMViterbiPackage/NotationHelper.cc)
add_library(edge_lib EMViterbiPackage/Edge.cc)
add_library(cypher_reader_lib CypherReader.cc)
add_library(graph_aid_lib EMViterbiPackage/GraphAid.cc)
add_library(trellis_aid_lib EMViterbiPackage/TrellisAid.cc)

add_executable(${PROJECT_NAME} Main.cc)
target_link_libraries(${PROJECT_NAME}
                      basic_helper_lib
                      notation_helper_lib
                      notation_lib
                      edge_lib
                      graph_aid_lib
                      cypher_reader_lib
                      trellis_aid_lib)

# Testing.
if (test)
#   add_executable(notation_test NotationTest.cc)
#   target_link_libraries(notation_test notation_lib)
endif()
