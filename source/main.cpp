// Include standard headers
#include <stdio.h>
#include <stdlib.h>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
using namespace glm;

#include <iostream>
#include <vector>
#include "common/shader.hpp"
#include "common/objloader.hpp"
#include "common/vboindexer.hpp"
#include "common/controls.hpp"
#include <glm/gtx/transform.hpp>
#include "View.hpp"
#include "Geometry.hpp"
#include "Fragment.hpp"
#include "Fragmenter.hpp"




int height=0;
int width=0;
float radius=0;
float densitysphere=0;
float densityline=0;
float maxpeak=0;
int numparts=0;
int lines=0;
int steps=0;
int seedLine=0;
int seedOrientation=0;
std::string ofile="";
int    m=6;           // initial vertex count; must be m>=3
float  niter=6;       // number of fractal iterations
float  jitter=0.3;      // relative jitter of midpoint along curve
float  amplitude=0.1;   // relative amplitude
float  decay=0.85;       // relative amplitude/jitter decay per iteration
float  smoothing=1;   // smoothing strength [0..1]

//Fragmenter::RandomFractureOptions  Fragmenter::defaultRandomFractureOptions = {    6, 6, 0.3, 0.3/3.0, 0.85, 1};

char* ofilepath;
int main(int argc, char *argv[] ){

    // std::cout <<  argc << std::endl;  //argv[0] is the program name
    srand(time(NULL));

    if(argc > 1)
    {
        for (int i = 1; i < argc; i++) {
            std::string option = argv[i];
            if (option == "-width") {
                //          std::cout << "here " << argv[i + 1]<< std::endl;
                width = atoi(argv[i + 1]);
            } else if (option == "-height") {
                //            std::cout << "height " << argv[i + 1]<< std::endl;
                
                height = atoi(argv[i + 1]);
            } else if (option == "-radius") {
                radius=atof(argv[i + 1]);
            }
            else if (option == "-densphere") {
                densitysphere=atof(argv[i + 1]);
            }
            else if (option == "-denline") {
                densityline=atof(argv[i + 1]);
            }
            else if (option == "-peak") {
                maxpeak=atof(argv[i + 1]);
            }
            else if (option == "-output") {
                ofile = argv[i + 1];
                ofilepath = argv[i+1];
            }
            else if (option == "-numparts") {
                numparts = atoi(argv[i + 1]);
            }
            else if (option == "-lines") {
                lines = atoi(argv[i + 1]);
            }
            else if (option == "-steps") {
                steps = atoi(argv[i + 1]);
            }
            else if (option == "-seedLine") {
                seedLine = atoi(argv[i + 1]);
            }
            else if (option == "-seedOrientation") {
                seedOrientation = atoi(argv[i + 1]);
            }
            else if (option == "-m") {
                m = atoi(argv[i + 1]);
            }
            else if (option == "-niter") {
                niter = atof(argv[i + 1]);
            }
            else if (option == "-jitter") {
                jitter = atof(argv[i + 1]);
            }
            else if (option == "-amplitude") {
                amplitude = atof(argv[i + 1]);
            }
            else if (option == "-decay") {
                decay = atof(argv[i + 1]);
            }
            else if (option == "-smoothing") {
                smoothing = atof(argv[i + 1]);
            }
        }
    }
    if (width==0) width = 500;
    if (height==0) height = 500;
    if (radius==0) radius = 5;
    if (densitysphere==0) densitysphere = 10;
    if (densityline==0) densityline = 20;
    if (maxpeak==0) maxpeak = 5;
    if (numparts==0) numparts = 10;
    if (lines==0) lines = 5;
    if (seedLine==0) seedLine = 10;
    if (seedOrientation==0) seedOrientation = 20;
    //initialise myview
    View* myview = new View(width, height);
    GLuint vertexarrayID = myview->initialiseVertexBuffer();

    //set user values to fragmenter
    Fragmenter::RandomFractureOptions useroptions{m, niter, jitter, amplitude, decay, smoothing};
    Fragmenter::setOptions(useroptions);
   
    //set user seeds
    Fragmenter::setSeedLine(seedLine);
    Fragment::setSeedOrientation(seedOrientation);

    Fragmenter* fragmenter = new Fragmenter(numparts,glm::vec3(1,0,0),radius,densityline,densitysphere,maxpeak,myview);

    fragmenter->fragment();
   // fragmenter->testCurve();

#if 1

    
    Geometry* geo = new Geometry(vertexarrayID,glm::vec3(0.7,0.7,0.7),GL_POINTS,GEO_SHAPE);
    geo->setDensity(densitysphere);
    geo->generateSphere();
    geo->calculateBoundingBox();
     geo->exportPoints();
    myview->addGeometry(geo);

    myview->calculateMaxBoundingBox();
    

    if (myview->initialise() != -1) {
        do {

        myview->update();
        //    jline->setNumSteps(steps);

        } // Check if the ESC key was pressed or the window was closed
        while (glfwGetKey(myview->window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
               glfwWindowShouldClose(myview->window) == 0);

    }    else {
        std::cout << "ERROR: The game has not been initialised correctly" << std::endl;
        return 0;
    }
#endif
}
