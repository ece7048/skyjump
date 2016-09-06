# CustomOptimization
#########################################################
# Source files

set(Target CustomOptimization)

file(GLOB INCS *.h pso/*.h)
file(GLOB SRCS *.cpp pso/*.cpp)
file(GLOB PSO_SRC pso/pso.cpp pso/testPSO.cpp)
file(GLOB PSO_INCS pSo/pso.h)
    
#########################################################
# Dependences

#########################################################
# Build library

addLibrary(
    KIT ${Target}
    LINKLIBS ${OPENSIMSIMBODY_LIBRARIES} 
    INCLUDES ${INCS}
    SOURCES ${SRCS}
)

add_executable(TestPSO ${PSO_SRC} ${PSO_INCS})
