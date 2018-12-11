#include "Geometry.hpp"



Geometry::Geometry(GLuint vertexarrayID, glm::vec3 colour, GLenum primitive, geo_type type): vertexArrayID(vertexArrayID), colour(colour), primitive(primitive), type(type){
	transformation = glm::mat4(1.0);
	bbtransformation = glm::mat4(1.0);
	visible = 1;
    numsteps = 1;

    vertexbuffer=0;
    colorbuffer=0;
    normalbuffer=0;
}
Geometry::~Geometry() {
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteBuffers(1, &colorbuffer);
}
void Geometry::setColor(glm::vec3 colort) {

	for (int i = 0; i < vertices.size(); i++) {
		colors.push_back(colort);
	}
}
void Geometry::storePointSet(std::vector<glm::vec3> verticest){

  /*  point_set.reserve (3); // For memory optimization
    for (std::vector<glm::vec3>::iterator it = vertices.begin() ; it != vertices.end(); ++it){
        glm::vec3 vertex = *it;
        
   //     std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
        point_set.insert(Point(double(vertex.x), double(vertex.y),double(vertex.z)));

    }*/
    
}

void Geometry::calculateBoundingBox() {
	boundingbox_min.x = boundingbox_max.x = vertices[0].x;
	boundingbox_min.y = boundingbox_max.y = vertices[0].y;
	boundingbox_min.z = boundingbox_max.z = vertices[0].z;
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].x < boundingbox_min.x) boundingbox_min.x = vertices[i].x;
		if (vertices[i].x > boundingbox_max.x) boundingbox_max.x = vertices[i].x;
		if (vertices[i].y < boundingbox_min.y) boundingbox_min.y = vertices[i].y;
		if (vertices[i].y > boundingbox_max.y) boundingbox_max.y = vertices[i].y;
		if (vertices[i].z < boundingbox_min.z) boundingbox_min.z = vertices[i].z;
		if (vertices[i].z > boundingbox_max.z) boundingbox_max.z = vertices[i].z;
	}

	//line1
	verticesBoundingBox.push_back(boundingbox_min); //back bottom left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_min.z)); //back bottom right
	//line2
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_min.z)); //back bottom right
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_max.y, boundingbox_min.z)); //back top right
	//line3
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_max.y, boundingbox_min.z)); //back top right
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_min.z)); //back top left
	//line4
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_min.z)); //back top left
	verticesBoundingBox.push_back(boundingbox_min); //back bottom left 
	//line 5
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_min.y, boundingbox_max.z)); //front bottom left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_max.z)); //front bottom right
	//line 6
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_max.z)); //front bottom right
	verticesBoundingBox.push_back(boundingbox_max); //front top right 
	//line 7
	verticesBoundingBox.push_back(boundingbox_max); //front top right 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_max.z)); //front top left 
	//line 8
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_max.z)); //front top left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_min.y, boundingbox_max.z)); //front bottom left 
	//line 9
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_min.y, boundingbox_min.z)); //front top left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_min.y, boundingbox_max.z)); //front bottom left 
	//line 10
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_min.z)); //front top left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_min.y, boundingbox_max.z)); //front bottom left 
	//line 11
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_min.z)); //front top left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_min.x, boundingbox_max.y, boundingbox_max.z)); //front bottom left 
	//line 12
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_max.y, boundingbox_min.z)); //front top left 
	verticesBoundingBox.push_back(glm::vec3(boundingbox_max.x, boundingbox_max.y, boundingbox_max.z)); //front bottom left 

}

std::vector<unsigned short> Geometry::getIndices() {
	return indices;
}
std::vector<glm::vec3> Geometry::getVertices() {
    
 //   std::cout << "$$$$$$$$$$$$$$$$" << vertices.size() << std::endl;
	return vertices;
}

std::vector<glm::vec3> Geometry::getVerticesReversed(){
    
    std::vector<glm::vec3> thevertices = vertices;
    std::reverse(thevertices.begin(),thevertices.end());
    return thevertices;
}
void Geometry::setVertices(std::vector<glm::vec3> verticest) {
    vertices=verticest;
}

std::vector<glm::vec3> Geometry::getColors() {
	//std::cout << colors.size() << std::endl;
	return colors;
}
std::vector<glm::vec2> Geometry::getUVs() {
	return uvs;
}
std::vector<glm::vec3> Geometry::getNormals() {
	return normals;
}
GLuint* Geometry::getVertexbuffer() {
	return &vertexbuffer;
}	
GLuint* Geometry::getColorbuffer() {
	return &colorbuffer;
}
GLuint* Geometry::getNormalbuffer() {
	return &normalbuffer;
}
GLuint* Geometry::getVertexBoundingBoxbuffer() {
	return &vertexbufferBoundingBox;
}
GLuint* Geometry::getColorBoundingBoxbuffer() {
	return &colorbufferBoundingBox;
}
std::vector<glm::vec3> Geometry::getColorsBoundingBox() {
	return colorsBoundingBox;
}
std::vector<glm::vec3> Geometry::getVerticesBoundingBox() {
	return verticesBoundingBox;
}


void Geometry::setInitialPosition(glm::vec3 initial_post) {
	initial_pos = initial_post;
	position = initial_pos;
	
}
void Geometry::updateBoundingBox() {
	boundingbox_min = boundingbox_min + position;
	boundingbox_max = boundingbox_max + position;
	bbtransformation = bbtransformation * glm::translate(position);
}
glm::vec3 Geometry::getInitialPosition() {
	return initial_pos;
}

glm::vec3 Geometry::getPosition() {
	return position;
}

void Geometry::setTransformation(glm::mat4 transformationt) {
	transformation = transformationt;
//	bbtransformation = bbtransformation * transformationt;

}
glm::mat4 Geometry::getTransformation() {
	return transformation;
}
glm::mat4 Geometry::getBoundingBoxTransformation() {
	return bbtransformation;
}
glm::vec3 Geometry::getBoundingBoxMin() {
	return boundingbox_min;
}
glm::vec3 Geometry::getBoundingBoxMax(){
	return boundingbox_max;

}


bool Geometry::isVisible() {
	return visible;
}
void Geometry::setVisible(bool is) {
	visible = is;
}

int Geometry::exportPoints(){
    

    /*
    // Writing result in OFF format
    int randomnum = rand() % 200;
    std::string filename = "points//points_"+std::to_string(randomnum)+".off";
   
    std::ofstream out(filename);
//    std::cout << "I am: " << point_set.size() << std::endl;
    if (!out || !CGAL::write_off_point_set (out, point_set))
    {
        return EXIT_FAILURE;
    }
    std::fstream file;
    file.open(filename);
    //file.seekp(0); // subtracts one from the buffer position ( now 18 )
   // file.write("", 1);
    file.close();
    return 0;
    return EXIT_SUCCESS;
 */
    
}
GLenum Geometry::getPrimitive(){
    return primitive;
}

void Geometry::generateSphere(){
    //  const double PI = 3.141592653589793238462643383279502884197;
    
    // int steps = density;
    float radius =1;
    for (float phi = 0; phi < glm::two_pi<float>(); phi+= glm::pi<float>()/density){// azimuth angle
        
        for (float theta = 0; theta < glm::pi<float>() ; theta+=glm::pi<float>()/density){//altitude/elevation angle
            double x = radius * cos(phi) * sin(theta);
            double y = radius * sin(phi) * sin(theta);
            double z = radius * cos(theta);
            vertices.push_back(glm::vec3(x,y,z));
            
        }
    }
    
    
    
    
    for (int i=0;i<vertices.size();i++){
        colors.push_back(colour);
    }
    for (int i=0;i<vertices.size();i++){
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
void Geometry::setDensity(float densityt){
    density = densityt;
}
geo_type Geometry::getType(){
    return type;
}
void Geometry::setType(geo_type typet){
    type = typet;
}
int Geometry::getNumSteps(){
    return numsteps;
}
void Geometry::setNumSteps(int numstepst){
    numsteps = numstepst;
}
