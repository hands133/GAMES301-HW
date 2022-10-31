# Introduction
Fork from GAMES301(https://ustc-gcl-f.github.io/code/index.html#sec_surface_framework)

Simplified CMakeLists.txt to use vcpkg to manage 3dparty for automatly building projects.


## GAMES 301 Homework

### Homework 2 : Analytic Eigensystems for Isotropic Distortion Energies (2D Symmetric Dirichlet Energy for now)

Prerequisites:
+ Tutte's Embedding structure

Target:
1. [**TODO**] Implement the basic structure of the eigensystem.

Usage:
 
```bash
git checkout hw2_EnergyEigensystems
./build_windows.bat
```
then run the subproject `SurfaceFrameworkCmake`.

Note

+ Relevant codes please refer to `Surface_Framework_Cmake/src/homeworks/EnergyEigensystems`:
  + `Util_EnergyEigensystems.h`,
  + `Util_EnergyEigensystems.cpp`.

Parameterization results

| Name | Format | Face Type | V | E | F | Boundaries | Storage |
| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|   cathead | OBJ | triangle |  131 |  378 |  248 | 1 |  8 KB |
|     Balls | OBJ | triangle |  547 | 1578 | 1032 | 1 | 26 KB |
| Bunnyhead | OBJ | triangle |  741 | 2188 | 1448 | 1 | 54 KB |
|      hand | OFF | triangle | 1558 | 4653 | 3096 | 1 | 97 KB |

### Gallery

**Under construction...**