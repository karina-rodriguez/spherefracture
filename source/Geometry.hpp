#pragma once
// Include GLM
#include <glm/glm.hpp>
#include <vector>
#include "common/objloader.hpp"
#include <glm/gtx/transform.hpp>
#include <string>     // std::string, std::to_string


#include <iostream>
#include <list>

//For pointset in CGAL
/*#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<unsigned char, 3> Color;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Property_map<Color> Color_map;
typedef Point_set::Property_map<FT> FT_map;*/
enum geo_type { GEO_SHAPE, GEO_FRAGMENT, GEO_INTERSECTION, GEO_PATH };

class Geometry {
    
    
protected:
    //attributes for drawing in OpenGL
    GLuint vertexArrayID;
	std::vector<unsigned short> indices;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> colors;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	GLuint vertexbuffer;
	GLuint colorbuffer;
	GLuint normalbuffer;

    glm::vec3 colour;
	glm::vec3 initial_pos;
	glm::vec3 position;
	glm::mat4 transformation;
	glm::mat4 bbtransformation;

    
	glm::vec3 boundingbox_min, boundingbox_max;
	bool visible;
    GLenum primitive;
    float density;
    geo_type type;
    int numsteps;

	std::vector<glm::vec3> verticesBoundingBox;
	std::vector<glm::vec3> colorsBoundingBox;
	GLuint vertexbufferBoundingBox;
	GLuint colorbufferBoundingBox;

    //attributes for storing information in CGAL
    //Point_set point_set;
   
    
public:
    int vertices_size;

    Geometry(GLuint vertexarrayID, glm::vec3 colour, GLenum primitive, geo_type type);
	~Geometry();

	std::vector<unsigned short> getIndices();
	
    //vertices
    std::vector<glm::vec3> getVertices();
    std::vector<glm::vec3> getVerticesReversed();
    void setVertices(std::vector<glm::vec3> verticest);

    
    std::vector<glm::vec3> getColors();
	std::vector<glm::vec2> getUVs();
	std::vector<glm::vec3> getNormals();
	GLuint* getVertexbuffer();
	GLuint* getNormalbuffer();
	GLuint* getColorbuffer();
	void setInitialPosition(glm::vec3 initial_pos);
	glm::vec3 getInitialPosition();
	glm::vec3 getPosition();
	void setColor(glm::vec3 colort);
	void calculateBoundingBox();
	void updateBoundingBox();
	glm::vec3 getBoundingBoxMin();
	glm::vec3 getBoundingBoxMax();
    GLenum getPrimitive();
    geo_type getType();
    void setType(geo_type typet);

	GLuint* getVertexBoundingBoxbuffer();
	GLuint* getColorBoundingBoxbuffer();
	std::vector<glm::vec3> getVerticesBoundingBox();
	std::vector<glm::vec3> getColorsBoundingBox();
	glm::mat4 getBoundingBoxTransformation();

	void setTransformation(glm::mat4 transformation);
	glm::mat4 getTransformation();
	bool isVisible();
	void setVisible(bool is);
    void generateSphere();
    void setDensity(float densityt);
    int getNumSteps();
    void setNumSteps(int numstepst);
    
    
    //functions of CGAL
    void storePointSet(std::vector<glm::vec3> verticest);
    int exportPoints();
    
        
        
        
};
