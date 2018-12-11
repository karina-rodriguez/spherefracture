#pragma once
// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>
#include <vector>

// Include GLM
#include <glm/glm.hpp>
using namespace glm;
#include <glm/gtx/transform.hpp>
#include "common/shader.hpp"
#include "Geometry.hpp"

class View {

private:
    
    
    GLuint programID;
    GLuint VertexArrayID;
    GLuint MatrixID;


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

    // Initial position : on +Z
    // Initial horizontal angle : toward -Z
    float horizontalAngle = 0.0f;
    // Initial vertical angle : none
    float verticalAngle = 0.0f;
    // Initial Field of View
    float initialFoV = 45.0f;
    
    float speed = 3.0f; // 3 units / second
    float mouseSpeed = 0.005f;
    glm::vec3 camposition = glm::vec3( 0, 0, 0);
    glm::vec3 direction  = glm::vec3( 0, 0, 0);
    
    
    
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
    void getMouseRotation();

    std::vector<Geometry*> getGeometries();
    void changeVisibility();

};

