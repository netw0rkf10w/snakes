project(application)
 
cmake_minimum_required(VERSION 3.8)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(ITK_DIR "/usr/local/itk/bin")
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
# use what QVTK built with
# SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
# SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
# SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "") 
#find_package(Qt5 REQUIRED)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

#include(${QT_USE_FILE})

find_package(Boost REQUIRED)

find_package(CGAL QUIET COMPONENTS Core REQUIRED)
include( ${CGAL_USE_FILE} )

# support for out-of-source build
INCLUDE_DIRECTORIES(
#  ${QT_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

 
set(SRC
    main.cpp
    application.h
    application.cpp
    view.h
    scene.h
    snake.h
    point.h
    dialog.h
    view.cpp
    scene.cpp
    snake.cpp
    point.cpp
    dialog.cpp
    stuffs.cpp
    stuffs.h
    application.ui
    dialog.ui
    Icons/icons.qrc)
 
#set(MOC application.h view.h scene.h snake.h point.h dialog.h stuffs.h)
#set(UI application.ui dialog.ui)
#set(RCC Icons/icons.qrc)
 
#qt5_wrap_cpp(MOC_GENERATED ${MOC})
#qt5_wrap_ui(UI_GENERATED ${UI})
#qt5_add_resources(RCC_GENERATED ${RCC})
 
#source_group("Form Files" FILES ${UI} )
#source_group("Resource Files" FILES ${RCC} ${EXE_ICON})
#source_group("Generated Files" FILES ${MOC_GENERATED}   ${UI_GENERATED}  ${RCC_GENERATED} )



add_executable(application
    ${SRC}
#    ${MOC_GENERATED}
#    ${UI_GENERATED}
#    ${RCC_GENERATED}
    )

#include_directories("/media/Khue/Dev. Software/C++libraries/boost")
 
target_link_libraries(application ${ITK_LIBRARIES}
#    ${QT_LIBRARIES}
    Qt5::Widgets
    Qt5::OpenGL

    )


####### To remove the dummy console window in MSVC, use the following configuration instead 
#add_executable(application WIN32
#    ${SRC}
#    ${MOC_GENERATED}   ${UI_GENERATED}  ${RCC_GENERATED} )
 
#target_link_libraries(application ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES} ${QT_QTMAIN_LIBRARY})

