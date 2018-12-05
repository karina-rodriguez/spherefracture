#include "View.hpp"

View::View(int width, int height):width(width),height(height){
    if (!glfwInit())
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        getchar();  
    }
    
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    // Open a window and create its OpenGL context
    window = glfwCreateWindow(width, height, "3D Viewer", NULL, NULL);
    if (window == NULL) {
        fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
        getchar();
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    // Hide the mouse and enable unlimited mouvement
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    
    // Set the mouse at the center of the screen
    glfwPollEvents();
    //we set the cursor at the centre so that it always start on our origin (0,0,0)
    glfwSetCursorPos(window, width / 2, height / 2);
    
    
    // Dark blue background
    glClearColor(1.0f, 1.0f, 1.0f, 0.5f);
    
    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    
    
    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
    }
    
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
    
    // Cull triangles which normal is not towards the camera
    //glEnable(GL_CULL_FACE);
    VertexArrayID = 1;
    counter = 0;
    FOV = 45;
    min = glm::vec3(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
    max = glm::vec3(0,0,0);
    glfwSetWindowPos(window, 0, 0);

}
View::~View() {
    // Cleanup VBO and shader
   
    glDeleteProgram(programID);
    glDeleteVertexArrays(1, &VertexArrayID);
    
    
    // Close OpenGL window and terminate GLFW
    glfwTerminate();
}


int View::initialise() {
    
   

    programID = LoadShaders("SimpleVertexShader.hlsl", "SimpleFragmentShader.hlsl");
    MatrixID = glGetUniformLocation(programID, "MVP");
   // ViewID = glGetUniformLocation(programID, "V");
    //ModelID = glGetUniformLocation(programID, "M");

    //LightID = glGetUniformLocation(programID, "LightPosition_worldspace");
    //CameraID = glGetUniformLocation(programID, "CameraPosition_worldspace");
    //lightPosition = glm::vec3(2, 2, 0);

    return 1;

}



void View::update() {

        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
        changeVisibility();

        glUseProgram(programID);
        glPolygonMode(GL_FRONT_AND_BACK, GL_TRIANGLES);
        for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
           Geometry* mygeo = *it;
            if (mygeo->isVisible()){

            float space=0.1;
          //  computeMatricesFromInputs(window, width, height,FOV,min.x-space,max.x+space,min.y-space,max.y+space);

            //computeMatricesFromInputs(<#GLFWwindow *window#>, <#int width#>, <#int height#>)
            //********Calculate the MVP matrix
            //***********PROJECTION*****************
           glm::mat4 Projection = glm::ortho(min.x-space,max.x+space,min.y-space,max.y+space,0.0f,200.0f);
                
           // glm::mat4 Projection = glm::ortho(min.x-space,max.x+space,min.y-space,max.y+space,0.0f,200.0f);

            //***********CAMERA*****************
            // Camera matrix
            
                
          //  std::cout << "position " << getPosition().x << ", " << getPosition().y << ", " << getPosition().z << std::endl;
            camposition = glm::vec3(boundingboxcentre.x, boundingboxcentre.y, distancez);
               // camposition = glm::vec3(0,0,0);

                getMouseRotation();
            glm::mat4 camera = glm::lookAt(
                                         //-glm::vec3(camposition.x, camposition.y, camposition.z), // Camera is at (0,0,-1), in World Space
                                           glm::vec3(camposition.x, camposition.y, camposition.z), // Camera is at (0,0,-1), in World Space
                                           
                                        glm::vec3(0,0,0), // and looks at the origin
                                         glm::vec3(0, 1, 0)
                                         );
               /* glm::mat4 camera1 = glm::lookAt(
                                               glm::vec3(0,0,0), // Camera is at (0,0,-1), in World Space
                                               direction, // and looks at the origin
                                               glm::vec3(0, 1, 0)
                                               );*/
                glm::mat4 View = camera; //* camera1; //* getViewMatrix();
            
            //***********MODEL*****************
            glm::mat4 Model = glm::mat4(1.0);//glm::rotate((float)90.0,glm::vec3(0.0,1.0,0.0));//glm::mat4(1.0);
            
            glm::mat4 mvp = Projection * View * Model;
            
            
            glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &mvp[0][0]);
            
            //********Add the Geometry
            // 1rst attribute buffer : vertices
            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, *mygeo->getVertexbuffer());
            glVertexAttribPointer(
                                  0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                                  3,                  // size
                                  GL_FLOAT,           // type
                                  GL_FALSE,           // normalized?
                                  0,                  // stride
                                  (void*)0            // array buffer offset
                                  );
            
            // 2nd attribute buffer : colors
            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, *mygeo->getColorbuffer());
            glVertexAttribPointer(
                                  1,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                                  3,                  // size
                                  GL_FLOAT,           // type
                                  GL_FALSE,           // normalized?
                                  0,                  // stride
                                  (void*)0            // array buffer offset
                                  );
            
            
            
            // 3rd attribute buffer : normals
            glEnableVertexAttribArray(2);
            glBindBuffer(GL_ARRAY_BUFFER, *mygeo->getNormalbuffer());
            glVertexAttribPointer(
                                  2,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
                                  3,                  // size
                                  GL_FLOAT,           // type
                                  GL_FALSE,           // normalized?
                                  0,                  // stride
                                  (void*)0            // array buffer offset
                                  );
            
            
            
            // Draw the geometry !
            glPointSize(5);

          //  glPolygonMode(GL_FRONT_AND_BACK, GL_TRIANGLES);
       ///     std::cout << mygeo->vertices_size << std::endl;
            glDrawArrays(mygeo->getPrimitive(), 0,mygeo->vertices_size); // 3 indices starting at 0 -> 1
          
            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(2);
            }
           

        }
        
        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();


}
void View::addGeometry(Geometry* mygeo){
    geometries.push_back(mygeo);
}
std::vector<Geometry*> View::getGeometries(){
    return geometries;
}
GLuint View::initialiseVertexBuffer(){
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);
    return VertexArrayID;
    
}
GLuint View::getVertexArrayID(){
    return VertexArrayID;
}
void View::calculateMaxBoundingBox(){
    distancez=0;
    glm::vec3 minbb(0,0,0);
    glm::vec3 maxbb(1000,1000,1000);
    for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
        Geometry* mygeo = *it;
        glm::vec3 ming = mygeo->getBoundingBoxMin();
        glm::vec3 maxg = mygeo->getBoundingBoxMax();
        
        if (ming.x<min.x) min.x = ming.x;
        if (ming.y<min.y) min.y = ming.y;
        if (ming.z<min.z) min.z = ming.z;
        
        if (maxg.x>max.x) max.x = maxg.x;
        if (maxg.y>max.y) max.y = maxg.y;
        if (maxg.z>max.z) max.z = maxg.z;
        
    }
    
    //std::cout << "min " << min.x << ", " << min.y << ", " << min.z << std::endl;
    //std::cout << "max " << max.x << ", " << max.y << ", " << max.z << std::endl;

    //calculate the ideal camera position
    //X Y are just in the centre of the camera
    boundingboxcentre.x = (min.x + max.x)/2;
    boundingboxcentre.y = (min.y + max.y)/2;
    boundingboxcentre.z = (min.z + max.z)/2;


    float distx = max.x-min.x;
    float disty = max.y-min.y;
    float diagonalg = sqrt(pow(distx,2)+pow(disty,2));
    
    // float distancez = diagonal/(atan(2*glm::radians(FOV/2)));
    //use alternative XY plane
    distancez = diagonalg/(atan(2*glm::radians(FOV/2)));
    std::cout << "distancez " << distancez << std::endl;

    
}
void View::changeVisibility(){

  //  std::cout << "START : " << geometries.size() << std::endl;

   
        //1 for fragments only
        if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS){
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
            //std::cout << "type:  " << mygeo->getType() << std::endl;

            if (mygeo->getType()==GEO_FRAGMENT) mygeo->setVisible(1);
            else mygeo->setVisible(0);
            }
        }
         //2 for paths only
        if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS){
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
            //    std::cout << "type:  " << mygeo->getType() << std::endl;
            if (mygeo->getType()==GEO_PATH) mygeo->setVisible(1);
            else mygeo->setVisible(0);
            }
        }
        //3 for intersections
        if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS){
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
             //   std::cout << "type:  " << mygeo->getType() << std::endl;

                if (mygeo->getType()==GEO_INTERSECTION) mygeo->setVisible(1);
                else mygeo->setVisible(0);
            }
        }

        // 0 for all
        if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS){
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
                mygeo->setVisible(1);
            }
        }
       // std::cout << mygeo->getNumSteps() << std::endl;

        // increase count
        // *************PAUSE*************
        static int oldState = GLFW_RELEASE;
        int newState = glfwGetKey(window, GLFW_KEY_I);
        if (newState == GLFW_RELEASE && oldState == GLFW_PRESS) {
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
                    mygeo->setNumSteps(mygeo->getNumSteps()+1);

            }
        }
        oldState = newState;
        
        static int oldState1 = GLFW_RELEASE;
        int newState1 = glfwGetKey(window, GLFW_KEY_D);
        if (newState1 == GLFW_RELEASE && oldState1 == GLFW_PRESS) {
            for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
                Geometry* mygeo = *it;
                std::cout << "type:  " << mygeo->getType() << std::endl;

            if (mygeo->getType()==GEO_FRAGMENT)
                if (mygeo->getNumSteps()>1) mygeo->setNumSteps(mygeo->getNumSteps()-1);
            }
        }
        oldState1 = newState1;
        
    
}
void View::getMouseRotation(){
    // glfwGetTime is called only once, the first time this function is called
    static double lastTime = glfwGetTime();
    
    // Compute time difference between current and last frame
    double currentTime = glfwGetTime();
    float deltaTime = float(currentTime - lastTime);
    
    
    // Get mouse position
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    
    // Reset mouse position for next frame
    glfwSetCursorPos(window, width / 2, height / 2);
    //std::cout << "****************** " << std::endl;
    
//    std::cout << "mouseSpeed: " <<  mouseSpeed << std::endl;
    
    // Compute new orientation
    horizontalAngle += mouseSpeed * float(width/2 - xpos );
    verticalAngle   += mouseSpeed * float(height/2 - ypos );
    
    //std::cout << "changes: " << (width / 2 - xpos) << ", " << (height / 2 - ypos) << std::endl;
    //std::cout << "horizontalAngle: " << horizontalAngle << std::endl;
    //std::cout << "verticalAngle: " << verticalAngle << std::endl;
    
    
    // Direction : Spherical coordinates to Cartesian coordinates conversion
    direction =glm::vec3(
                         cos(verticalAngle) * sin(horizontalAngle),
                         sin(verticalAngle),
                         cos(verticalAngle) * cos(horizontalAngle)
                         );
    //    std::cout << "direction: " << direction.x << ", " << direction.y << ", " <<
    //    direction.z << std::endl;
    
    // Right vector
    glm::vec3 right = glm::vec3(
                                sin(horizontalAngle - 3.14f/2.0f),
                                0,
                                cos(horizontalAngle - 3.14f/2.0f)
                                );
    //std::cout << "right: " << right.x << ", " << right.y << ", " <<
    //    right.z << std::endl;
    
    // Up vector
    glm::vec3 up = glm::cross( right, direction );
    //std::cout << "position inital: " << position.x << ", " << position.y << ", " << position.z << std::endl;
    
    // Move forward
    if (glfwGetKey( window, GLFW_KEY_UP ) == GLFW_PRESS){
        camposition += direction * deltaTime * speed;
    }
    // Move backward
    if (glfwGetKey( window, GLFW_KEY_DOWN ) == GLFW_PRESS){
        camposition -= direction * deltaTime * speed;
    }
    // Strafe right
    if (glfwGetKey( window, GLFW_KEY_RIGHT ) == GLFW_PRESS){
        camposition += right * deltaTime * speed;
    }
    // Strafe left
    if (glfwGetKey( window, GLFW_KEY_LEFT ) == GLFW_PRESS){
        camposition -= right * deltaTime * speed;
    }
  //  std::cout << "position after changes: " << camposition.x << ", " << camposition.y << ", " << camposition.z << std::endl;
    
//std::cout << "position after changes: " << direction.x << ", " << direction.y << ", " << direction.z << std::endl;
    
        
    // For the next frame, the "last time" will be "now"
    lastTime = currentTime;
    
}
