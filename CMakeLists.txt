# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( breakSphere )


#cmake_minimum_required(VERSION 2.8.11)

cmake_minimum_required(VERSION 3.1.0)   # for CMAKE_CXX_STANDARD
set (CMAKE_CXX_STANDARD 11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

find_package(OpenGL REQUIRED)


if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()


# Compile external dependencies 
add_subdirectory (external)

# On Visual 2005 and above, this module can set the debug working directory
cmake_policy(SET CMP0026 OLD)


include_directories(
	external/glfw-3.1.2/include/GLFW/
	external/glm-0.9.7.1/
	external/glew-1.13.0/include/
    external/assimp-3.0.1270/include/
	.
)


set(ALL_LIBS
	${OPENGL_LIBRARY}
	glfw
	GLEW_1130
)

add_definitions(
#	-DTW_STATIC
#	-DTW_NO_LIB_PRAGMA
#	-DTW_NO_DIRECT3D
	-DGLEW_STATIC
	-D_CRT_SECURE_NO_WARNINGS
)

#FILE(GLOB SRC_FILES "source/*.hlsl")

set_source_files_properties(${SRC_FILES} PROPERTIES VS_TOOL_OVERRIDE "None")

# breakSphere
add_executable(breakSphere 
	${SRC_FILES}		
	source/main.cpp
	source/View.cpp
	source/View.hpp
	source/Geometry.cpp
	source/Geometry.hpp
	source/Fragment.cpp
	source/Fragment.hpp
	source/Fragmenter.cpp
	source/Fragmenter.hpp
	source/SimpleFragmentShader.hlsl
	source/SimpleVertexShader.hlsl
	common/objloader.cpp
    common/objloader.hpp
    common/shader.cpp
	common/shader.hpp
)

# working directories
set_target_properties(breakSphere PROPERTIES COMPILE_DEFINITIONS "USE_ASSIMP")
target_link_libraries(breakSphere
	${ALL_LIBS}
	${CGAL_LIBRARIES}
	${CGAL_3RD_PARTY_LIBRARIES}
	assimp
)




set_target_properties(breakSphere PROPERTIES XCODE_ATTRIBUTE_CONFIGURATION_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source/")

#message(${CMAKE_GENERATOR})
if(NOT ${CMAKE_GENERATOR} STREQUAL "Xcode")
add_custom_command(
   TARGET breakSphere POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/breakSphere${CMAKE_EXECUTABLE_SUFFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/source/"
    COMMENT "Copying breakSphere files..."

)
endif()



SOURCE_GROUP(common REGULAR_EXPRESSION ".*/common/.*" )
SOURCE_GROUP(shaders REGULAR_EXPRESSION ".*/.*hlsl*" )



