# Framework for surface mesh processing （C++Framework for  Assignments of GAMES301 ）

**The program , which runs in the Qt application, is used to load meshes and render using OpenGL.**


## Dependent External Libraries
* [Qt](https://www.qt.io/), Recommended version: 5.13.0
## Alternative External Libraries
* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/), Recommended version: the latest 8.1(at Oct. 2020)

## Usage
* Use vcpkg to install qt
* Set the VCPKG_ROOT environment variable to the vcpkg root directory
* Use cmake to build 
```
cd Surface_Framework_Cmake
mkdir build 
cmake -S . -B build
```

Open **SurfaceFrameworkCmake.sln**, select **SurfaceFrameworkCmake** in build as launch project, and run.

## Supported File Formats

.obj .off .ply .stl