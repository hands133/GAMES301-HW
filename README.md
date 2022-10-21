# Introduction
Fork from GAMES301(https://ustc-gcl-f.github.io/code/index.html#sec_surface_framework)

Simplified CMakeLists.txt to use vcpkg to manage 3dparty for automatly building projects.


## GAMES 301 Homework

### Homework 1 : Tutte's Embedding (Planar Parameterization)

Target:
1. [**Done**] Implement the basic structure of the Tutte's Embedding method.
2. [**Done**] Calculate the weight of lambda with Tutte's (uniform) and Floater's (shape-preserving) respectively.
3. [**Done**] Visualize the UV mapping.

Usage:
 
```bash
git checkout hw1
./build_windows.bat
```
then run the subproject `SurfaceFrameworkCmake`.

Note

+ I've add tabs of Tutte's weight and Floater Weight, but the UV mapping is only visible under `smooth rendering` mode.
+ If you want to visualize the 2D parameterization results, uncomment the line 683 in `MeshViewerWidget.cpp`.
+ If you want to change the 2D parameterization boundary shape, please refer to line 631 to line 635 in `MeshViewerWidget.cpp`.
+ If you want to add customized boundary shape, please add enums in `enum class UVBoundaryType` in line 376, and implement the corresponding codes calculating UV in switch-case branch in line 393 in `MeshViewerWidget.cpp`.
+ If you've got good ideas, raise Issue and create Pull requests please.