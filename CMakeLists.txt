cmake_minimum_required(VERSION 3.2)

set (CMAKE_CXX_STANDARD 17)
SET(ITK_DIR "C:\\Libs\\InsightToolkit-5.2.1\\build64")
set(Boost_Root "C:\\Libs\\boost_1_72_0")
set(BOOST_DIR "C:\\Libs\\boost_1_72_0")

set(CMAKE_TOOLCHAIN_FILE "C:/Libs/vcpkg/scripts/buildsystems/vcpkg.cmake"
  CACHE STRING "Vcpkg toolchain file")

project(CT2Surface)
 
 
string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
string(REPLACE "/MD" "/MTd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
 
set (CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Default build type: Release" FORCE)
set(CMAKE_MFC_FLAG 1) 

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

find_package(ITK)
include(${ITK_USE_FILE})

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

find_package(Boost COMPONENTS program_options)

find_package(CGAL)

#set(sources
#   CT2Surface.cpp
#    )

#set(headers
#   Misc.h
#    )
#add_compile_options(${sources} ${headers} /W4 /WX)

add_executable(CT2Surface CT2Surface.cpp Misc.h)


set_property(TARGET CT2Surface PROPERTY
  MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Release>:Release>")

TARGET_LINK_LIBRARIES( CT2Surface LINK_PUBLIC ${Boost_LIBRARIES} ${ITK_LIBRARIES} CGAL::CGAL)
