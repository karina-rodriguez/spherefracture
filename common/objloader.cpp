#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>

#include <glm/glm.hpp>

#include "objloader.hpp"
#include <assimp/IOSystem.hpp>
#include <assimp/Exporter.hpp>



// Very, VERY simple OBJ loader.
// Here is a short list of features a real function would provide : 
// - Binary files. Reading a model should be just a few memcpy's away, not parsing a file at runtime. In short : OBJ is not very great.
// - Animations & bones (includes bones weights)
// - Multiple UVs
// - All attributes should be optional, not "forced"
// - More stable. Change a line in the OBJ file and it crashes.
// - More secure. Change another line and you can inject code.
// - Loading from memory, stream, etc

bool loadOBJ(
	const char * path, 
	std::vector<glm::vec3> & out_vertices, 
	std::vector<glm::vec2> & out_uvs,
	std::vector<glm::vec3> & out_normals
){
	printf("Loading OBJ file %s...\n", path);

	std::vector<unsigned int> vertexIndices, uvIndices, normalIndices;
	std::vector<glm::vec3> temp_vertices; 
	std::vector<glm::vec2> temp_uvs;
	std::vector<glm::vec3> temp_normals;


	FILE * file = fopen(path, "r");
	if( file == NULL ){
		printf("Impossible to open the file ! Are you in the right path ? See Tutorial 1 for details\n");
		getchar();
		return false;
	}

	while( 1 ){

		char lineHeader[128];
		// read the first word of the line
		int res = fscanf(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File. Quit the loop.

		// else : parse lineHeader
		
		if ( strcmp( lineHeader, "v" ) == 0 ){
			glm::vec3 vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
			temp_vertices.push_back(vertex);
		}else if ( strcmp( lineHeader, "vt" ) == 0 ){
			glm::vec2 uv;
			fscanf(file, "%f %f\n", &uv.x, &uv.y );
			uv.y = -uv.y; // Invert V coordinate since we will only use DDS texture, which are inverted. Remove if you want to use TGA or BMP loaders.
			temp_uvs.push_back(uv);
		}else if ( strcmp( lineHeader, "vn" ) == 0 ){
			glm::vec3 normal;
			fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
			temp_normals.push_back(normal);
		}else if ( strcmp( lineHeader, "f" ) == 0 ){
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
			//int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
			int matches = fscanf(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0], &vertexIndex[1], &normalIndex[1], &vertexIndex[2], &normalIndex[2]);
			if (matches != 6){
				printf("File can't be read by our simple parser :-( Try exporting with other options\n");
				fclose(file);
				return false;
			}
			vertexIndices.push_back(vertexIndex[0]);
			vertexIndices.push_back(vertexIndex[1]);
			vertexIndices.push_back(vertexIndex[2]);
			//uvIndices    .push_back(uvIndex[0]);
			//uvIndices    .push_back(uvIndex[1]);
			//uvIndices    .push_back(uvIndex[2]);
			normalIndices.push_back(normalIndex[0]);
			normalIndices.push_back(normalIndex[1]);
			normalIndices.push_back(normalIndex[2]);
		}else{
			// Probably a comment, eat up the rest of the line
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}
	std::cout << vertexIndices.size() << std::endl;
	// For each vertex of each triangle
	for( unsigned int i=0; i<vertexIndices.size(); i++ ){

		// Get the indices of its attributes
		unsigned int vertexIndex = vertexIndices[i];
		//unsigned int uvIndex = uvIndices[i];
		unsigned int normalIndex = normalIndices[i];
		
		// Get the attributes thanks to the index
		glm::vec3 vertex = temp_vertices[ vertexIndex-1 ];
		//glm::vec2 uv = temp_uvs[ uvIndex-1 ];
		glm::vec3 normal = temp_normals[ normalIndex-1 ];
		
		// Put the attributes in buffers
		std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
		out_vertices.push_back(vertex);
		//out_uvs     .push_back(uv);
		out_normals .push_back(normal);
	
	}
	fclose(file);
	return true;
}


#ifdef USE_ASSIMP // don't use this #define, it's only for me (it AssImp fails to compile on your machine, at least all the other tutorials still work)

// Include AssImp
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

bool loadAssImpMultipleMeshes(const char * path, std::vector<struct MyMesh> & myMeshes) {
	Assimp::Importer importer;
	std::cout << "1 "  << std::endl;

	const aiScene* scene = importer.ReadFile(path, 0/*aiProcess_JoinIdenticalVertices | aiProcess_SortByPType*/);
	std::cout << "2 " << std::endl;

	//std::cout << scene->mMeshes[2];
	if (!scene) {
		fprintf(stderr, importer.GetErrorString());
		getchar();
		return false;
	}
	std::cout << "3 " << scene->mNumMeshes << std::endl;

	// For each mesh
	for (unsigned int n = 0; n < scene->mNumMeshes; ++n)
	{

		std::cout << "4 " << std::endl;

		struct MyMesh* aMesh = new MyMesh;
		const aiMesh* mesh = scene->mMeshes[n]; // In this simple example code we always use the 1rst mesh (in OBJ files there is often only one anyway)
												// create array with faces
												// have to convert from Assimp format to array
		
			/*unsigned int *faceArray;
		faceArray = (unsigned int *)malloc(sizeof(unsigned int) * mesh->mNumFaces * 3);
		unsigned int faceIndex = 0;
		std::cout << "5 " << std::endl;

		for (unsigned int t = 0; t < mesh->mNumFaces; ++t) {
			const aiFace* face = &mesh->mFaces[t];
			
			memcpy(&faceArray[faceIndex], face->mIndices, 3 * sizeof(unsigned int));
			faceIndex += 3;
		}
		aMesh->numFaces = mesh->mNumFaces;*/
		std::cout << "5 " << mesh->mNumVertices << std::endl;
		
		aMesh->numFaces = mesh->mNumFaces;

		// Fill vertices positions
		aMesh->vertices.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
			aiVector3D pos = mesh->mVertices[i];
			aMesh->vertices.push_back(glm::vec3(pos.x, pos.y, pos.z));
		}

		
		// Fill vertices texture coordinates
		/*uvs.reserve(mesh->mNumVertices);
		for(unsigned int i=0; i<mesh->mNumVertices; i++){
		aiVector3D UVW = mesh->mTextureCoords[0][i]; // Assume only 1 set of UV coords; AssImp supports 8 UV sets.
		uvs.push_back(glm::vec2(UVW.x, UVW.y));
		}*/

		// Fill vertices normals
		aMesh->normals.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
			aiVector3D n = mesh->mNormals[i];
			aMesh->normals.push_back(glm::vec3(n.x, n.y, n.z));
		}


		// Fill face indices
		aMesh->indices.reserve(3 * mesh->mNumFaces);
		for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
			// Assume the model has only triangles.
			aMesh->indices.push_back(mesh->mFaces[i].mIndices[0]);
			aMesh->indices.push_back(mesh->mFaces[i].mIndices[1]);
			aMesh->indices.push_back(mesh->mFaces[i].mIndices[2]);
		}
		std::cout << "indices " << aMesh->vertices.size() << std::endl;
		myMeshes.push_back(*aMesh);

	}
	return true;
}
bool loadAssImpMeshes(
	const char * path,
	std::vector<unsigned short> & indices,
	std::vector<glm::vec3> & vertices,
	std::vector<glm::vec2> & uvs,
	std::vector<glm::vec3> & normals,
	std::vector<unsigned int> & list_vertices

) {

	Assimp::Importer importer;

	const aiScene* scene = importer.ReadFile(path, 0/*aiProcess_JoinIdenticalVertices | aiProcess_SortByPType*/);

	//std::cout << scene->mMeshes[2];
	if (!scene) {
		fprintf(stderr, importer.GetErrorString());
		getchar();
		return false;
	}

	list_vertices.reserve(scene->mNumMeshes);

	// For each mesh
	for (unsigned int n = 0; n < scene->mNumMeshes; ++n)
	{

		const aiMesh* mesh = scene->mMeshes[n]; // In this simple example code we always use the 1rst mesh (in OBJ files there is often only one anyway)
		list_vertices.push_back(mesh->mNumVertices);
												// Fill vertices positions
		vertices.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
			aiVector3D pos = mesh->mVertices[i];
			vertices.push_back(glm::vec3(pos.x, pos.y, pos.z));
		}

		// Fill vertices texture coordinates
		/*uvs.reserve(mesh->mNumVertices);
		for(unsigned int i=0; i<mesh->mNumVertices; i++){
		aiVector3D UVW = mesh->mTextureCoords[0][i]; // Assume only 1 set of UV coords; AssImp supports 8 UV sets.
		uvs.push_back(glm::vec2(UVW.x, UVW.y));
		}*/

		// Fill vertices normals
		normals.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
			aiVector3D n = mesh->mNormals[i];
			normals.push_back(glm::vec3(n.x, n.y, n.z));
		}


		// Fill face indices
		indices.reserve(3 * mesh->mNumFaces);
		for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
			// Assume the model has only triangles.
			indices.push_back(mesh->mFaces[i].mIndices[0]);
			indices.push_back(mesh->mFaces[i].mIndices[1]);
			indices.push_back(mesh->mFaces[i].mIndices[2]);
		}
	}

	// The "scene" pointer will be deleted automatically by "importer"
	return true;
}


bool loadAssImp(
	const char * path,
	std::vector<unsigned short> & indices,
	std::vector<glm::vec3> & vertices,
	std::vector<glm::vec2> & uvs,
	std::vector<glm::vec3> & normals
	
) {

	Assimp::Importer importer;

	const aiScene* scene = importer.ReadFile(path, 0/*aiProcess_JoinIdenticalVertices | aiProcess_SortByPType*/);

	//std::cout << scene->mMeshes[2];
	if (!scene) {
		fprintf(stderr, importer.GetErrorString());
		getchar();
		return false;
	}
	

		const aiMesh* mesh = scene->mMeshes[0]; // In this simple example code we always use the 1rst mesh (in OBJ files there is often only one anyway)


												// Fill vertices positions
		vertices.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i<mesh->mNumVertices; i++) {
			aiVector3D pos = mesh->mVertices[i];
			vertices.push_back(glm::vec3(pos.x, pos.y, pos.z));
		}

		// Fill vertices texture coordinates
		/*uvs.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i<mesh->mNumVertices; i++) {
			aiVector3D UVW = mesh->mTextureCoords[0][i]; // Assume only 1 set of UV coords; AssImp supports 8 UV sets.
			uvs.push_back(glm::vec2(UVW.x, UVW.y));
		}*/

		// Fill vertices normals
		normals.reserve(mesh->mNumVertices);
		for (unsigned int i = 0; i<mesh->mNumVertices; i++) {
			aiVector3D n = mesh->mNormals[i];
			normals.push_back(glm::vec3(n.x, n.y, n.z));
		}


		// Fill face indices
		indices.reserve(3 * mesh->mNumFaces);
		for (unsigned int i = 0; i<mesh->mNumFaces; i++) {
			// Assume the model has only triangles.
			indices.push_back(mesh->mFaces[i].mIndices[0]);
			indices.push_back(mesh->mFaces[i].mIndices[1]);
			indices.push_back(mesh->mFaces[i].mIndices[2]);
		}


	// The "scene" pointer will be deleted automatically by "importer"
	return true;
}

bool exportScene( const std::string& pFile,
                 std::vector<glm::vec3> thevertices,
                 const std::vector<glm::vec3> thenormals,
                 const std::vector<int> theindices)
{
    aiScene scene;
    scene.mRootNode = new aiNode();
    
    //material
    scene.mMaterials = new aiMaterial*[ 1 ];
    scene.mMaterials[ 0 ] = nullptr;
    scene.mNumMaterials = 1;
    scene.mMaterials[ 0 ] = new aiMaterial();
    
    
    scene.mMeshes = new aiMesh*[ 1 ];
    scene.mMeshes[ 0 ] = nullptr;
    scene.mNumMeshes = 1;
    
    scene.mMeshes[ 0 ] = new aiMesh();
    scene.mMeshes[ 0 ]->mMaterialIndex = 0;
    
    scene.mRootNode->mMeshes = new unsigned int[ 1 ];
    scene.mRootNode->mMeshes[ 0 ] = 0;
    scene.mRootNode->mNumMeshes = 1;
    
    auto pMesh = scene.mMeshes[ 0 ];
    
    std::cout << " VERTICES SIZE " << thevertices.size() << " " << thenormals.size() << " " << theindices.size() << std::endl;
    
    
    /*
    
    //Default Fill Location Vector
    int draw_order[36] =
    {
        0,2,1,      2,3,1,
        1,3,5,      3,7,5,
        5,7,4,      7,6,4,
        4,6,0,      6,2,0,
        4,0,5,      0,1,5,
        2,6,3,      6,7,3
    };
    
    glm::vec3 data[8] =
    {
        glm::vec3(-1.0f/2.0f,1.0f/2.0f,1.0f/2.0f),
        glm::vec3(1.0f/2.0f,1.0f/2.0f,1.0f/2.0f),
        glm::vec3(-1.0f/2.0f,-1.0f/2.0f,1.0f/2.0f),
        glm::vec3(1.0f/2.0f,-1.0f/2.0f,1.0f/2.0f),
        glm::vec3(-1.0f/2.0f,1.0f/2.0f,-1.0f/2.0f),
        glm::vec3(1.0f/2.0f,1.0f/2.0f,-1.0f/2.0f),
        glm::vec3(-1.0f/2.0f,-1.0f/2.0f,-1.0f/2.0f),
        glm::vec3(1.0f/2.0f,-1.0f/2.0f,-1.0f/2.0f)
    };
    
    for(int i = 0; i < 36; i++)
    {
        vertices.push_back(data[draw_order[i]]);
    }
    
    //Default Fill Normal Vector
    for(int i = 0; i < 36; i++)
    {
        if(i < 6)       {normals.push_back(glm::vec3(0,0,1));}
        else if(i < 12) {normals.push_back(glm::vec3(1,0,0));}
        else if(i < 18) {normals.push_back(glm::vec3(0,0,-1));}
        else if(i < 24) {normals.push_back(glm::vec3(-1,0,0));}
        else if(i < 30) {normals.push_back(glm::vec3(0,1,0));}
        else if(i < 36) {normals.push_back(glm::vec3(0,- 1,0));}
    }
    
    //Default Fill UV Vector
    for(int i = 0; i < 6; i++)
    {
        uvs.push_back(glm::vec2(0,1));
        uvs.push_back(glm::vec2(0,0));
        uvs.push_back(glm::vec2(1,1));
        uvs.push_back(glm::vec2(0,0));
        uvs.push_back(glm::vec2(1,0));
        uvs.push_back(glm::vec2(1,1));
    }*/
    
    //const auto& vVertices = thevertices;
    
    pMesh->mVertices = new aiVector3D[ thevertices.size() ];
    pMesh->mNormals = new aiVector3D[ thevertices.size() ];
    
    pMesh->mNumVertices = thevertices.size();
    
 //   pMesh->mTextureCoords[ 0 ] = new aiVector3D[ vVertices.size() ];
 //   pMesh->mNumUVComponents[ 0 ] = vVertices.size();
    
    int j = 0;
    for ( auto itr = thevertices.begin(); itr != thevertices.end(); ++itr )
    {
        pMesh->mVertices[ itr - thevertices.begin() ] = aiVector3D( thevertices[j].x, thevertices[j].y, thevertices[j].z );
        pMesh->mNormals[ itr - thevertices.begin() ] = aiVector3D( thenormals[j].x, thenormals[j].y, thenormals[j].z );
  //      pMesh->mTextureCoords[0][ itr - vVertices.begin() ] = aiVector3D( uvs[j].x, uvs[j].y, 0 );
        j++;
    }
    
    pMesh->mFaces = new aiFace[ thevertices.size() / 3 ];
    pMesh->mNumFaces = (unsigned int)(thevertices.size() / 3);
    
    int k = 0;
    for(int i = 0; i < (thevertices.size() / 3); i++)
    {
        aiFace &face = pMesh->mFaces[i];
        face.mIndices = new unsigned int[3];
        face.mNumIndices = 3;
        
        face.mIndices[0] = k;
        face.mIndices[1] = k+1;
        face.mIndices[2] = k+2;
        k = k + 3;
    }
    
    //mExportFormatDesc->id is "collada"  and mFilePath is "C:/Users/kevin/Desktop/myColladaFile.dae"
   // mAiExporter->Export(scene, mExportFormatDesc->id, mFilePath);
    
    Assimp::Exporter* exp;
    exp = new Assimp::Exporter();
    const char* id = "stl";
    exp->Export(&scene, id, pFile);
 
    return 1;
}
#endif
