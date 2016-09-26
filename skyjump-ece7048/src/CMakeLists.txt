CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
set(TARGET SkyJumpOptimization)

file(GLOB INCS *.h)
file(GLOB SRCS *.cpp)

add_executable(${TARGET} ${SRCS} ${INCS})
target_link_libraries(${TARGET} ${OPENSIMSIMBODY_LIBRARIES})
