#ifndef CONTROLS_HPP
#define CONTROLS_HPP

void computeMatricesFromInputs(GLFWwindow* window, int width, int height, float FOV, float leftw, float rightw, float bottomw, float topw);
glm::vec3 getPosition();
glm::mat4 getViewMatrix();
glm::mat4 getProjectionMatrix();

#endif
