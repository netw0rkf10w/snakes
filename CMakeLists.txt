project(application)
 
cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_BUILD_TYPE Debug)
 
set(ITK_DIR "/usr/local/itk/bin")
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
# use what QVTK built with
# SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
# SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
# SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "") 
find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})
 
find_package(CGAL QUIET COMPONENTS Core REQUIRED)
include( ${CGAL_USE_FILE} )

# support for out-of-source build
INCLUDE_DIRECTORIES(
  ${QT_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

 
set(SRC main.cpp application.h application.cpp view.h scene.h snake.h point.h dialog.h view.cpp scene.cpp snake.cpp point.cpp dialog.cpp stuffs.cpp stuffs.h)
 
set(MOC application.h view.h scene.h snake.h point.h dialog.h stuffs.h)
set(UI application.ui dialog.ui)
set(RCC Icons/icons.qrc)
 
qt4_wrap_cpp(MOC_GENERATED ${MOC})
qt4_wrap_ui(UI_GENERATED ${UI})
qt4_add_resources(RCC_GENERATED ${RCC})
 
source_group("Form Files" FILES ${UI} )
source_group("Resource Files" FILES ${RCC} ${EXE_ICON})
source_group("Generated Files" FILES ${MOC_GENERATED}   ${UI_GENERATED}  ${RCC_GENERATED} )



add_executable(application
    ${SRC}
    ${MOC_GENERATED}   ${UI_GENERATED}  ${RCC_GENERATED} )

include_directories("/media/Khue/Dev. Software/C++libraries/boost")
 
target_link_libraries(application ${ITK_LIBRARIES} ${QT_LIBRARIES})


####### To remove the dummy console window in MSVC, use the following configuration instead 
#add_executable(application WIN32
#    ${SRC}
#    ${MOC_GENERATED}   ${UI_GENERATED}  ${RCC_GENERATED} )
 
#target_link_libraries(application ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES} ${QT_QTMAIN_LIBRARY})

