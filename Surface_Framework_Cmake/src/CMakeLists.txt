cmake_minimum_required(VERSION 3.18)

find_package(Qt5 REQUIRED 
        COMPONENTS Widgets Core Gui OpenGL)

set(EXTRA_LINK_LIB
    Qt5::Widgets 
    Qt5::Core
    Qt5::Gui
    Qt5::OpenGL
    opengl32
    glu32  
)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp *.h *.ui)
file(GLOB_RECURSE resoures CONFIGURE_DEPENDS *.qrc)
qt5_add_resources(source ${resoures})

file (GLOB GUI_FILE
    main.cpp
    MeshViewer/MainViewerWidget.*
    MeshViewer/MeshParamWidget.*
    surfacemeshprocessing.h
    surfacemeshprocessing.cpp
)

file (GLOB MeshViewer_FILE
    MeshViewer/InteractiveViewerWidget.*
    MeshViewer/MeshViewerWidget.*
    MeshViewer/QGLViewerWidget.*
)

file (GLOB PolyMesh_FILE
    PolyMesh/IOManger.cpp
    PolyMesh/PolyMesh.cpp
    PolyMesh/PolyMesh_Base.cpp
    PolyMesh/PolyMeshIterators.cpp
    include/Math/Algorithms.h
    include/Math/Massage.h
    include/Math/MPoint3.h
    include/Math/Mrgb.h
    include/Math/MVector3.h
    include/PolyMesh/IOManger.h
    include/PolyMesh/MemoryPool.h
    include/PolyMesh/PolyMesh.h
    include/PolyMesh/PolyMesh_Base.h
    include/PolyMesh/PolyMeshIterators.h
)

source_group(GUI FILES ${GUI_FILE})
source_group(MeshViewer FILES ${MeshViewer_FILE})
source_group(PolyMesh FILES ${PolyMesh_FILE})

add_executable(${PROJECT_NAME}
    ${source}
)

target_link_libraries(${PROJECT_NAME} PRIVATE ${EXTRA_LINK_LIB})
