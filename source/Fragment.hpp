#pragma once
// Include GLM
#include <glm/glm.hpp>
#include "Geometry.hpp"

#include <iostream>
#include <list>


class Fragment : public Geometry {

    

public:
   
    Fragment(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, std::vector<glm::vec3> verticest);
    ~Fragment();


};
