#include "Fragment.hpp"

template <typename T>
inline T const& max (T const& a, T const& b, T const& c) {
    
    double maxv = ( a < b ) ? b : a;
    return ( ( maxv < c ) ? c : maxv );
 }
const double  Fragment::epsilon = 1e-6;

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
//! Create a plane defined by a point and a normal.
/*
 *
 *  Returns true if function succed, and populate the centroid and normal
 *  of the plane
 *  Solution taken from: https://www.ilikebigbits.com/2015_03_04_plane_from_points.html
 *
 */
int Fragment::createPlane(const std::vector<glm::vec3> vertices){
   
   // std::vector<glm::vec3> vertices = getVertices();

    //at least 3 points required
    if  (vertices.size()< 3)
        return 0;
    glm::vec3 sum(0.0, 0.0, 0.0);
    
    for (int i=0;i<vertices.size();i++) {
        sum += vertices[i];
    }
    
    std::cout << "sum: " << sum.x << ", " << sum.y << ", " << sum.z << std::endl;
  
    glm::vec3 centroid = glm::vec3(sum.x/ vertices.size(), sum.y/ vertices.size(), sum.z/ vertices.size());
    std::cout << "centroid: " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;
    
    
    // Calc full 3x3 covariance matrix, excluding symmetries:
    double xx = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yy = 0.0;
    double yz = 0.0;
    double zz = 0.0;
   
    for (int i=0;i<vertices.size();i++) {
        glm::vec3 r = vertices[i] - centroid;
        xx += r.x * r.x;
        xy += r.x * r.y;
        xz += r.x * r.z;
        yy += r.y * r.y;
        yz += r.y * r.z;
        zz += r.z * r.z;
    }
    
    double det_x = yy*zz - yz*yz;
    double det_y = xx*zz - xz*xz;
    double det_z = xx*yy - xy*xy;
    
    double det_max = max(det_x, det_y, det_z);
    std::cout << "Determinants: " << det_x << ", " << det_y << ", " << det_z << ", " << det_max << std::endl;
    
    
    if (det_max <= 0.0) {
        return 0; // The points don't span a plane
    }
    
    // Pick path with best conditioning:
    glm::vec3 dir;
    if (det_max == det_x) dir = glm::vec3(det_x, xz*yz - xy*zz,xy*yz - xz*yy);
    else{
        if (det_max == det_y) dir = glm::vec3(xz*yz - xy*zz, det_y, xy*xz - yz*xx);
        else {
            dir = glm::vec3( xy*yz - xz*yy,xy*xz - yz*xx, det_z);
        }
    }
    glm::vec3 norm1 = glm::normalize(dir);

    //theplane = {norm1, centroid};
    return 1;
}

//! Intersection between a line containing the point in the fragment and a plane
/*! Plane is defined by the centroid and normal
 *
 *  Returns true if there is an interesection point and return the point itslef
 *  If the lines are parallel it returns false, assuming we do not need to use
 *  that point.
 *  It also makes the assumption that the line is made by the point in the fragment
 *  and the centre of the sphere which is (0,0,0)
*/
int Fragment::checkIntersectionwithPlane(glm::vec3 point, glm::vec3& result){
    double lambda;
 
    //the initial point is (0,0,0) as we take the centre of the sphere for the vector
    //which defines the line for intersection
    static const glm::vec3 p0(0,0,0);
    //check that the lines are not parallel to each other
    if ((glm::dot(point,theplane.normal)>epsilon)
        ||(glm::dot(point,theplane.normal)<epsilon)){
        lambda = glm::dot((theplane.centroid - point),theplane.normal)/
                 glm::dot((point-p0),theplane.normal);
        
        std::cout << "intersect? " << glm::dot(point,theplane.normal) <<  " lambda: "  << lambda << std::endl;
       

        std::cout << "resulting point: " << (point.x*lambda)+point.x << ", " << (point.y*lambda)+point.y << ", " << (point.z*lambda)+point.z  << std::endl;

    }else {
        return 0;
        
    }
    result = glm::vec3((point.x*lambda)+point.x,(point.y*lambda)+point.y,(point.z*lambda)+point.z);
   // result = point;
    return 1;
    
}

