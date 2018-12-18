#ifndef OBJLOADER_H
#define OBJLOADER_H

// TO READ VARIOUS MESHES
// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <glfw3.h>

// Include GLM
#include <glm/glm.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>



struct MyMesh {

	int numFaces;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> faces;
	std::vector<glm::vec2>  uvs;
	std::vector<glm::vec3>  normals;

	std::vector<unsigned short> indices;

};

bool loadAssImpMultipleMeshes(const char * path, std::vector<struct MyMesh> & myMeshes);
bool loadAssImpMeshes(
	const char * path,
	std::vector<unsigned short> & indices,
	std::vector<glm::vec3> & vertices,
	std::vector<glm::vec2> & uvs,
	std::vector<glm::vec3> & normals,
	std::vector<unsigned int> & list_vertices

);

bool loadOBJ(
	const char * path, 
	std::vector<glm::vec3> & out_vertices, 
	std::vector<glm::vec2> & out_uvs, 
	std::vector<glm::vec3> & out_normals
);


bool loadAssImp(
	const char * path, 
	std::vector<unsigned short> & indices,
	std::vector<glm::vec3> & vertices,
	std::vector<glm::vec2> & uvs,
	std::vector<glm::vec3> & normals
);

bool exportScene( const std::string& pFile,
                 const std::vector<glm::vec3> thevertices,
                 const std::vector<glm::vec3> thenormals,
                 const std::vector<int> theindices);

bool myOwnExportSceneSTL( const std::string& pFile,
                         std::vector<glm::dvec3> thevertices,
                         const std::vector<glm::dvec3> thenormals);



#endif
