set(TARGET SkyJumpOptimization)

add_subdirectory(customoptimization)
include_directories(customoptimization)

file(GLOB INCS *.h)
file(GLOB SRCS *.cpp)

add_executable(${TARGET} ${SRCS} ${INCS})
target_link_libraries(${TARGET} ${OPENSIMSIMBODY_LIBRARIES} CustomOptimization)
