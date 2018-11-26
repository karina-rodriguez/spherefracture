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
    max = glm::vec3(std::numeric_limits<float>::min(),std::numeric_limits<float>::min(),std::numeric_limits<float>::min());
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

            /////////***************************************
            /////////************FIRST VIEWPORT***********
            /////////***************************************
            glViewport(0, 0, width*2, height*2);

            float space=0.1;
            computeMatricesFromInputs(window, width, height,FOV,min.x-space,max.x+space,min.y-space,max.y+space);

            //computeMatricesFromInputs(<#GLFWwindow *window#>, <#int width#>, <#int height#>)
            //********Calculate the MVP matrix
            //***********PROJECTION*****************
            glm::mat4 Projection = getProjectionMatrix();
         //   glm::perspective(glm::radians(FOV), 1.0f / 1.0f, 0.1f, 100.0f);
            //getProjectionMatrix();
            //***********CAMERA*****************
            // Camera matrix
            
            glm::vec3 cameraPosition(boundingboxcentre.x, boundingboxcentre.y, distancez*15);
            
            glm::mat4 camera = glm::lookAt(
                                         glm::vec3(cameraPosition.x, cameraPosition.y, cameraPosition.z), // Camera is at (0,0,-1), in World Space
                                         //glm::vec3(0, 0, 2),
                                         glm::vec3(boundingboxcentre.x, boundingboxcentre.y, boundingboxcentre.z), // and looks at the origin
                                         glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
                                         );
            glm::mat4 View = camera * getViewMatrix();
            //***********MODEL*****************
           // glm::mat4 Model = glm::rotate((float)counter,glm::vec3(0.0,1.0,0.0));
            //counter+=0.01;
            //glm::scale(glm::vec3(0.05,0.05,0.05));
            glm::mat4 Model = glm::mat4(1.0);//glm::rotate((float)90.0,glm::vec3(0.0,1.0,0.0));//glm::mat4(1.0);
            
            glm::mat4 mvp = Projection * View * Model;
            
            
            glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &mvp[0][0]);
      //      glUniformMatrix4fv(ModelID, 1, GL_FALSE, &View[0][0]);
       //     glUniformMatrix4fv(ViewID, 1, GL_FALSE, &View[0][0]);

         //   glUniform3f(LightID, lightPosition.x, lightPosition.y, lightPosition.z);
//            glUniform3f(CameraID, cameraPosition.x, cameraPosition.y, cameraPosition.z);

        

            
            
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
    for (std::vector<Geometry*>::iterator it = geometries.begin() ; it != geometries.end(); ++it){
        Geometry* mygeo = *it;
        glm::vec3 ming = mygeo->getBoundingBoxMin();
        glm::vec3 maxg = mygeo->getBoundingBoxMax();
        
        if (ming.x<min.x) min.x = ming.x;
        if (ming.y<min.y) min.y = ming.y;
        if (maxg.x>max.x) max.x = maxg.x;
        if (maxg.y>max.y) max.y = maxg.y;
        //calculate the ideal camera position
        glm::vec3 boundingboxcentreg;
        boundingboxcentreg.x = (ming.x + maxg.x)/2;
        float distx = maxg.x-ming.x;
        boundingboxcentreg.y = (ming.y + maxg.y)/2;
        float disty = maxg.y-ming.y;
        float distz = maxg.z-ming.z;
        
        boundingboxcentreg.z = (ming.z + maxg.z)/2;
        boundingboxcentre = boundingboxcentreg;
        float diagonalg = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
        
        // float distancez = diagonal/(atan(2*glm::radians(FOV/2)));
        //use alternative XY plane
        distancez = diagonalg/(atan(2*glm::radians(FOV/2)));
        //glm::mat4 Projection = glm::ortho(-10.0f,10.0f,-10.0f,10.0f,0.0f,100.0f); // In world coordinates
        //left right bottom top
        
    }
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
