cmake_minimum_required (VERSION 3.0)
project(sergey-dalitz-fit)

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

include(cmake/settings.cmake)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install bin to:  " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# require a fairly recent ROOT version
# this is needed by all subdirectories
find_package(ROOT "5.30" REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})
# including them as SYSTEM prevents
# many false-positive warnings
include_directories(SYSTEM ${ROOT_INCLUDE_DIR})

