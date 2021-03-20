cmake_minimum_required(VERSION 3.16)
project(nDFVM)

 # c++11 std
 set(CMAKE_CXX_STANDARD 11)
 # We need a CMAKE_DIR with some code to find external dependencies
 set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
# set Sources

 file(GLOB SOURCES "src/*.cpp" "src/*.h")

 add_executable(nDFVM ${SOURCES})
