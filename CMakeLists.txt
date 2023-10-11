cmake_minimum_required(VERSION 3.23)
project(cecher)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL REQUIRED)

add_executable(cecher cecher.cpp)

target_compile_definitions(cecher PRIVATE -DNDEBUG)

target_link_libraries(cecher CGAL::CGAL)

option(USE_COEFFICIENTS "enable support for coefficients in a prime field Z/pZ (default p=2)" OFF)
option(EXPLICIT_CENTERS "enable variant of the minimal enclosing sphere algorithm for low ambient dimension" OFF)

if(USE_COEFFICIENTS)
	target_compile_definitions(cecher PRIVATE USE_COEFFICIENTS)
endif()

if(EXPLICIT_CENTERS)
	target_compile_definitions(cecher PRIVATE EXPLICIT_CENTERS)
endif()