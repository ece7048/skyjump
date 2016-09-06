CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

project(SkyJump)
#########################################################
# project

set(PROJECT_LIBRARIES "" CACHE INTERNAL "Contains all project libraries")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(cmake/CMakeHelpers.cmake)

# path
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} CACHE PATH "")
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR} CACHE PATH "")

# Set the variables that would have otherwise been set by installer
set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_INCLUDEDIR include)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_DOCDIR doc)

# build type
set(CMAKE_BUILD_TYPE Release CACHE string "Debug or Release build" FORCE)

# build options
option(USE_UNITY_BUILD "For fast building" on)

# for grouping into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Settings
set(BASE_DIR "${PROJECT_SOURCE_DIR}")
set(INI_PATH "${PROJECT_SOURCE_DIR}/data/setup.ini")
set(CMAKE_EXTRA_INCLUDE_FILES src/Settings.h)
configure_file("${PROJECT_SOURCE_DIR}/src/Settings.h.in" "${PROJECT_SOURCE_DIR}/src/Settings.h")
set(CMAKE_EXTRA_INCLUDE_FILES)

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.

if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(PLATFORM_ABI x64)
else()
    set(PLATFORM_ABI x86)
endif()

set(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" 
    CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

# VS warnings
if(MSVC)
    add_definitions("/W3 /D_CRT_SECURE_NO_WARNINGS /wd4996 /wd4800 /wd4267 /wd4244 /wd4251")
endif()

#################################################################################################
# Find OpenSim

find_package(OpenSim REQUIRED)
if(OpenSim_FOUND)
    message(STATUS "OpenSim found at ${OPENSIM_ROOT_DIR}")
    include_directories(${OPENSIMSIMBODY_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Can't find OpenSim at ${OPENSIM_ROOT_DIR}")
endif()

#################################################################################################
# subdirectories

add_subdirectory(src)



