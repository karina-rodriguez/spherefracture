#include "Fragment.hpp"

Fragment::Fragment(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, std::vector<glm::vec3> verticest): Geometry(vertexarrayIDT,colour,primitive, type) {
 
    
   /* vertices.push_back(glm::vec3(-0.5, 0, 0));
    vertices.push_back(glm::vec3(0.5, 0, 0));
    vertices.push_back(glm::vec3(0, 1, 0));
    */
    vertices = verticest;
    for (int i=0;i<vertices.size();i++){
        colors.push_back(colour);
    }
    for (int i=0;i<vertices.size();i++){
        //        normals.push_back(glm::vec3(0.0, 0.0, 10.0));
        normals.push_back(vertices[i]-glm::vec3(0,0,0));
        
    }
    
    
    storePointSet(vertices);
    
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);
    vertices_size =vertices.size();
    
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), &colors[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
    
}

Fragment::~Fragment() {
}
void Fragment::setLevel(int levelt){
    level = levelt;
}
int Fragment::getLevel(){
    return level;
}
