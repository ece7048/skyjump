# Sky Jump Optimization

Status: Untested

Find a set of control points to perform the best possible high jump. The following code is based on [SimTK](http://simtk-confluence.stanford.edu:8080/display/OpenSim/Sky+High%3A+Coordinating+Muscles+for+Optimal+Jump+Performance) and is adjusted to perform an optimization based on a set of control nodes and an objective function similar to [1].

[1] Anderson F.C., el al., A Dynamic Optimization Solution for Vertical Jumping in Three Dimensions, Computer Methods in Biomechanics and Biomedical Engineering, 2:3, 201-32, 1999

# Dependency

It is based on OpenSim 3.3. The user must define the OPENSIM_HOME variable and to add OPENSIM_HOME/bin to path. 

# Building

Run CMake in the root directory.

# Running

You can configure your simulation through data/setup.ini, where the path configured by CMake through src/Setting.h.