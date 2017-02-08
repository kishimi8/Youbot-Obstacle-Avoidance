# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 2.8.3)
project( assign2 )

find_package(catkin REQUIRED COMPONENTS
CGAL QUIET COMPONENTS 
  roscpp
  turtlesim
  geometry_msgs
)
 
# CGAL and its components
if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

catkin_package()
# include helper file
include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories(
  ${catkin_INCLUDE_DIRS}
)


add_executable(planner src/planner.cpp)
# add_executable(planner src/polygon_test.cpp)

target_link_libraries(planner
   ${catkin_LIBRARIES}
 )
