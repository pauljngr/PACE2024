cmake_minimum_required(VERSION 3.5)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)  
set (CMAKE_CXX_STANDARD 11)

project(OSCM2024 LANGUAGES CXX)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".dylib")
set(CMAKE_BUILD_TYPE "Release")

find_package(CBC REQUIRED)
find_package(CLP REQUIRED)
find_package(COINUTILS REQUIRED)
find_package(OSI REQUIRED)
find_package(CGL REQUIRED)

add_executable(oscm oscm.cpp CbcInterface.cpp)

target_include_directories(oscm PRIVATE ${CBC_INCLUDE_DIR} ${CLP_INCLUDE_DIR} ${COINUTILS_INCLUDE_DIR} ${OSI_INCLUDE_DIR} ${CGL_INCLUDE_DIR})

target_link_libraries(oscm PRIVATE ${CBC_LIBRARIES} ${CLP_LIBRARIES} ${COINUTILS_LIBRARIES} ${OSI_LIBRARIES} ${CGL_LIBRARIES})

target_compile_options(oscm PRIVATE -Wall -O3 -DNDEBUG) #-fsanitize=address
target_link_options(oscm PRIVATE -Wall -O3 -DNDEBUG) # -DNDEBUG
