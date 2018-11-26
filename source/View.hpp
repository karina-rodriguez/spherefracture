#pragma once
// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>
#include <vector>

// Include GLM
#include <glm/glm.hpp>
using namespace glm;
#include "common/objloader.hpp"
#include "common/vboindexer.hpp"
#include <glm/gtx/transform.hpp>
#include "common/shader.hpp"
#include "common/controls.hpp"
#include "common/text2D.hpp"
#include "Geometry.hpp"

class View {

private:
    
    
    GLuint programID;
    GLuint VertexArrayID;
    GLuint MatrixID;
    GLuint ModelID;
    GLuint ViewID;
    GLuint LightID;
    GLuint CameraID;

    //screen size
    int width;
    int height;
     float counter;
    std::vector<Geometry*> geometries;
    glm::vec3 lightPosition;
    
    float FOV;
    glm::vec3 boundingboxcentre;
    glm::vec3 min;
    glm::vec3 max;
    float distancez;

    
public:
    GLFWwindow* window;
    View(int width, int height);
    ~View();
    int initialise();
    void calculateMaxBoundingBox();
    GLuint initialiseVertexBuffer();
    GLuint getVertexArrayID();
    void update();
    void addGeometry(Geometry* mygeo);

    std::vector<Geometry*> getGeometries();
    void changeVisibility();

};

