cmake_minimum_required(VERSION 3.10)

# Set the project name
project(Inverse)

# Tell CMake to use CTest extension
enable_testing()

# Include the subdirectory containing headers
include_directories(src/include)
include_directories(thirdparty/muparser/include)
include_directories(thirdparty/cnpy)

# Add the cnpy libray
add_library(cnpy SHARED "thirdparty/cnpy/cnpy.cpp")

# Default mode is Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Demand lots of warnings and all warnings as errors
# set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Werror")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") 
# add_compile_options(-Wall -Wextra -pedantic -Werror)


# Add source files needed for my target "transfer"
add_executable(transfer 
    src/main.cpp 
    src/exporter.cpp 
    src/solver.cpp 
    src/mesh.cpp 
    src/config.cpp
    thirdparty/muparser/src/muParser.cpp
    thirdparty/muparser/src/muParserBase.cpp
    thirdparty/muparser/src/muParserBytecode.cpp
    thirdparty/muparser/src/muParserCallback.cpp
    thirdparty/muparser/src/muParserError.cpp
    thirdparty/muparser/src/muParserTokenReader.cpp)

# Link to OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(transfer OpenMP::OpenMP_CXX)
endif()

# Link the cnpy library to my target
target_link_libraries(transfer cnpy)

# Find and link the ZLIB to my target (needed by cnpy)
find_package(ZLIB)
target_link_libraries(transfer ZLIB::ZLIB)

# A simple test that works perfectly for sds
add_test(NAME run_sds COMMAND transfer ../src/config/test.cfg)
