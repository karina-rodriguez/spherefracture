#pragma once
// Include GLM
#include <glm/glm.hpp>
#include "Geometry.hpp"

#include <iostream>
#include <list>


class JaggedLine : public Geometry {

    int randomness;
    float density;
    float maxpeak;

    void  initFromVertices(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type);
public:
   
    JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, const std::vector<glm::vec3> &poly);
    JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, float density, float maxpeak, bool anticlockwise);
    ~JaggedLine();
    JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, float density,float maxpeak,bool anticlockwise, glm::vec3 rot);
};
