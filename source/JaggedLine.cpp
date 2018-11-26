#include "JaggedLine.hpp"

JaggedLine::JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, float density,float maxpeak,bool anticlockwise): Geometry(vertexarrayIDT,colour, primitive, type), density(density), maxpeak(maxpeak) {
//    const double PI = 3.141592653589793238462643383279502884197;
    //generate a random rotation by giving a random vector
    float radius = 1;
    float randomx = ((float)rand()/RAND_MAX * 2) - 1;
    float randomy = ((float)rand()/RAND_MAX * 2) - 1;
    float randomz = ((float)rand()/RAND_MAX * 2) - 1;
    std::cout << "random : " << randomx << ", " << randomy <<", " << randomz << std::endl;

    glm::mat4 rotation = glm::rotate((float)90.0,glm::vec3(randomx,randomy,randomz));
  //  glm::mat4 rotation = glm::rotate((float)90,glm::vec3(1,0,0));

    //phi is angle in XZ plane, theta is angle in XY plane
    // iterate through all aximuth angles, but only for one altitude: PI/2
    float  phi_step =  glm::pi<float>()/density;
    //check if we need to do one phi_step less to avoid having two points similar to each other
    for (float phi = 0; phi < glm::two_pi<float>()-phi_step; phi+=phi_step){// azimuth angle
        float theta = glm::pi<float>()/2;//altitude angle
        
        float randomvalue = ((float)rand()/RAND_MAX * (2*(0.6*phi_step))) - (0.6*phi_step);
        float ph = phi + randomvalue;
        
        /*
        double x;
        double xtmp = radius * cos(ph) * sin(theta);
        if (anticlockwise) x = -xtmp;
        else*/
        double x = radius * cos(ph) * sin(theta);
        double y = radius * sin(ph) * sin(theta);
        double ztmp = radius * cos(theta);
        //get a random value between -maxpeak*2 and maxpeak*2
        float randomzvalue = ((float)rand()/RAND_MAX * (maxpeak*2)) - maxpeak;
        double z = ztmp + randomzvalue;
       //  double z = radius * cos(theta);
//std::cout << "here : " << x << ", " << y <<", " << z << std::endl;
        //-x makes it go anticlockwise
        glm::vec4 p(-x,y,z,0);
   //     std::cout << "p : " << p.x << ", " << p.y <<", " << p.z << std::endl;
        glm::vec4 newp = p * rotation;
       std::cout << "newp : " << newp.x << ", " << newp.y <<", " << newp.z << std::endl;
        vertices.push_back(glm::vec3(newp.x,newp.y,newp.z));
//  vertices.push_back(glm::vec3(p.x,p.y,p.z));
    }
    
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
    
    std::cout << "**********$$$$$$$$$$$$$$$$" << vertices.size() << std::endl;

    
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), &colors[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
}
JaggedLine::~JaggedLine() {
}
//testing jagged line
JaggedLine::JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, float density,float maxpeak,bool anticlockwise, glm::vec3 rot): Geometry(vertexarrayIDT,colour, primitive, type), density(density), maxpeak(maxpeak) {
    //    const double PI = 3.141592653589793238462643383279502884197;
    //generate a random rotation by giving a random vector
    float radius = 1;
    float randomx = ((float)rand()/RAND_MAX * 2) - 1;
    float randomy = ((float)rand()/RAND_MAX * 2) - 1;
    float randomz = ((float)rand()/RAND_MAX * 2) - 1;
    std::cout << "random : " << randomx << ", " << randomy <<", " << randomz << std::endl;
    int num=0;
    glm::mat4 rotation = glm::rotate((float)90.0,rot);
    //  glm::mat4 rotation = glm::rotate((float)90,glm::vec3(1.0,0.0,0.0));
    
    //phi is angle in XZ plane, theta is angle in XY plane
    // iterate through all aximuth angles, but only for one altitude: PI/2
    float  phi_step =  glm::pi<float>()/density;
    //check if we need to do one phi_step less to avoid having two points similar to each other
    for (float phi = 0; phi < glm::two_pi<float>()-phi_step; phi+=phi_step){// azimuth angle
        float theta = glm::pi<float>()/2;//altitude angle
        
        float randomvalue = ((float)rand()/RAND_MAX * (2*(0.6*phi_step))) - (0.6*phi_step);
        float ph = phi;// + randomvalue;
        
        /*
         double x;
         double xtmp = radius * cos(ph) * sin(theta);
         if (anticlockwise) x = -xtmp;
         else*/
        double x = radius * cos(ph) * sin(theta);
        double y = radius * sin(ph) * sin(theta);
        double ztmp = radius * cos(theta);
        //get a random value between -maxpeak*2 and maxpeak*2
        float randomzvalue = ((float)rand()/RAND_MAX * (maxpeak*2)) - maxpeak;
        double z = ztmp;// + randomzvalue;
        //  double z = radius * cos(theta);
        //std::cout << "here : " << x << ", " << y <<", " << z << std::endl;
        //-x makes it go anticlockwise
        glm::vec4 p(-x,y,z,0);
        //     std::cout << "p : " << p.x << ", " << p.y <<", " << p.z << std::endl;
        glm::vec4 newp = p * rotation;
        std::cout << "newp : " << newp.x << ", " << newp.y <<", " << newp.z << std::endl;
        vertices.push_back(glm::vec3(newp.x,newp.y,newp.z));
        //  vertices.push_back(glm::vec3(p.x,p.y,p.z));
       // if (num==5) break;
     //   num++;
    }
    
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
    
    std::cout << "**********$$$$$$$$$$$$$$$$" << vertices.size() << std::endl;
    
    
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), &colors[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
}
