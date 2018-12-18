#include "Fragment.hpp"

// [TW:] The following definition returned references into the stack,
// which is why I replaced it by an alternative definition...
/* template <typename T>
inline T const& max (T const& a, T const& b, T const& c) {
    double maxv = ( a < b ) ? b : a;
    return ( ( maxv < c ) ? c : maxv );
    } */
// ... and here is that alternative definition:
template <typename T>
inline T const& max (T const& a, T const& b, T const& c) { return std::max(std::max(a, b), c); }

const double  Fragment::epsilon = 1e-6;

Fragment::Fragment(GLuint vertexarrayIDT,
                   glm::vec3 colour,
                   GLenum primitive,
                   geo_type type,
                   std::vector<glm::dvec3> verticest)
    : Geometry(vertexarrayIDT,colour,primitive, type)
{
 
    
   /* vertices.push_back({-0.5, 0, 0});
    vertices.push_back({0.5, 0, 0});
    vertices.push_back({0, 1, 0});
    */
    vertices = std::vector<glm::vec3>(verticest.begin(), verticest.end());
    colors = std::vector<glm::vec3>(vertices.size(), colour);
    for (int i=0;i<vertices.size();i++){
        //        normals.push_back({0.0, 0.0, 10.0});
        normals.push_back(vertices[i]-glm::vec3(0,0,0));
        
    }
    
    maxdistancecentroidfragment=0;
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
 *  Returns true if function succeeds, and populate the centroid and normal
 *  of the plane
 *  Code based on: https://www.ilikebigbits.com/2015_03_04_plane_from_points.html
 *  Main modifications include the generation of 3 planes: one which fits the vertices, a
 *  second one which is beyond the bounding box and one close to the centre.
 *
 *
 */
int Fragment::createPlanes(const std::vector<glm::dvec3> vertices, const double close, const double far){
   
   // std::vector<glm::dvec3> vertices = getVertices();

    //at least 3 points required
    if  (vertices.size()< 3)
        return 0;
    glm::dvec3 sum(0.0, 0.0, 0.0);
    
    for (int i=0;i<vertices.size();i++) {
        sum += vertices[i];
    }
    
   // std::cout << "sum: " << sum.x << ", " << sum.y << ", " << sum.z << std::endl;
  
    glm::dvec3 centroid = (1.0/vertices.size()) * sum;
    
    
    // Calc full 3x3 covariance matrix, excluding symmetries:
    double xx = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yy = 0.0;
    double yz = 0.0;
    double zz = 0.0;
   
    for (int i=0;i<vertices.size();i++) {
        glm::dvec3 r = vertices[i] - centroid;
        double distance = glm::distance(vertices[i],centroid);
        
      //  std::cout << "distance: " << distance << std::endl;
        if (distance>maxdistancecentroidfragment) {
            maxdistancecentroidfragment=distance;
            maxdistancecentroidfragmentvertex =vertices[i];
        }
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
  //  std::cout << "Determinants: " << det_x << ", " << det_y << ", " << det_z << ", " << det_max << std::endl;
    
    
    if (det_max <= 0.0) {
        return 0; // The points don't span a plane
    }
    
    // Pick path with best conditioning:
    glm::dvec3 dir;
    if (det_max == det_x) dir = glm::dvec3(det_x, xz*yz - xy*zz,xy*yz - xz*yy);
    else{
        if (det_max == det_y) dir = glm::dvec3(xz*yz - xy*zz, det_y, xy*xz - yz*xx);
        else {
            dir = glm::dvec3( xy*yz - xz*yy,xy*xz - yz*xx, det_z);
        }
    }
    

    glm::dvec3 norm = glm::normalize(dir);
    glm::dvec3 normcentroid = glm::normalize(centroid);

    //std::cout << "centroid: " << centroid.x << ", " << centroid.y << ", " << centroid.z << std::endl;
    //std::cout << "normcentroid: " << normcentroid.x << ", " << normcentroid.y << ", " << normcentroid.z << std::endl;
   // std::cout << "normcentroid * 1.5: " << normcentroid.x*1.5 << ", " << normcentroid.y*1.5 << ", " << normcentroid.z*1.5 << std::endl;




    //create a plane with the normal and centroid
    plane = {norm, centroid};
    glm::dvec3 centroidclosesttocentre(normcentroid.x*close,
                    normcentroid.y*close,normcentroid.z*close);
    closestplane = {norm, centroidclosesttocentre};
    
    glm::dvec3 centroidfurthesttocentre(normcentroid.x*far,
                    normcentroid.y*far,normcentroid.z*far);
    furthestplane = {norm, centroidfurthesttocentre};

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
int Fragment::checkIntersectionwithPlane(glm::dvec3 point, Plane plane, glm::dvec3& result)
{
    double lambda;
 
    //the initial point is (0,0,0) as we take the centre of the sphere for the vector
    //which defines the line for intersection
    static const glm::dvec3 p0(0,0,0);
    //check that the lines are not parallel to each other
    if ((glm::dot(point,plane.normal)>epsilon)
        ||(glm::dot(point,plane.normal)<epsilon)){
        lambda = glm::dot((plane.centroid - point),plane.normal)/
                 glm::dot((point-p0),plane.normal);
        
     //   std::cout << "intersect? " << glm::dot(point,theplane.normal) <<  " lambda: "  << lambda << std::endl;
       

    //    std::cout << "resulting point: " << (point.x*lambda)+point.x << ", " << (point.y*lambda)+point.y << ", " << (point.z*lambda)+point.z  << std::endl;

    }else {
        return 0;
        
    }
    result = lambda * point + point;
   // result = point;
    return 1;
    
}

Plane Fragment::getPlane(){
    return plane;
}
Plane Fragment::getClosestPlane(){
    return closestplane;
}
Plane Fragment::getFurthestPlane(){
    return furthestplane;
}
//! Returns a transformation for bringing the plane to XY plane
/*
 *
 * Simply computes the rotation required for facing the XY plane
 *
 */
glm::mat4 Fragment::getTransformationForPolygoninXYPlane(Plane theplane){
    
    glm::vec3 normal = theplane.normal;
    glm::vec3 centroid = theplane.centroid;

    glm::vec3 targetaxis = {0,0,1};
    //find the axis of rotation through the cross-product of the two normals
    glm::vec3 axisofrotation = glm::cross(targetaxis,glm::normalize(normal));
    
    std::cout << "axisofrotation: " << axisofrotation.x << ", " << axisofrotation.y << ", " << axisofrotation.z << std::endl;

    double tetha = acos(glm::dot(glm::normalize(normal),targetaxis));
    std::cout << "TETHA: " << tetha << std::endl;
    std::cout << "TETHA: " << (float)glm::degrees(tetha) << std::endl;
    //rotate by an angle in the axis of rotation so that normal aligns with (0,0,1)
    glm::mat4 rot = glm::rotate((float)tetha, glm::normalize(axisofrotation));
    
    return rot;
    
    
}

//! Function to create an STL file for the fragment
/*
 *
 * Performs several steps:
 *    - Create three planes: one which fits the vertices, one further beyond the bounding box and
 *  one close to the bounding box. The values for these are closestvalue and furtherstvalue
 *    - Computes the intersection between each vertex in the fragment point and the plane
 *  to create a polygon which lies on the plane.
 *    - check all points are actually on the plane. Although this does not affect anything else for the
 *  moment.
 *    - Transforms the plane to the XY plane by requesting the transformation and applying it to each point
 *    - Computes a polygon partition using CGAL - this requires of simple polygons!
 *
 */
void Fragment::createSTL(int counterfile){
    static const double closestvalue = 0.001;
    static const double furtherstvalue = 4;

    //view->addGeometry(tmpqueue.front());
    //create the planes for intersection, one near the centre, the other far away from the radius of the sphere
    createPlanes(getVerticesD(), closestvalue, furtherstvalue);

    
   /* for (int n=0;n<getVertices().size();n++){
        
        std::cout << getVertices()[n].x << " " << getVertices()[n].y << " " << getVertices()[n].z << std::endl;
        glm::vec3 pointcloseplane,pointfarplane;
        
        //get the points of the side closest to the centre of the sphere
        checkIntersectionwithPlane(getVertices()[n],getClosestPlane(),pointcloseplane);
        closeplanepoints.push_back(pointcloseplane);
        
        //get the points of the side furthest to the centre of the sphere
        checkIntersectionwithPlane(getVertices()[n],getFurthestPlane(),pointfarplane);
        farplanepoints.push_back(pointfarplane);
        
        
    }*/
    

    
    //check if they are actually on the plane
  /*  bool allpointsonplane=1;
    for (int n=0;n<farplanepoints.size();n++){
        
        glm::vec3 v1(farplanepoints[n]-getFurthestPlane().centroid);
        if(!(glm::dot(v1,getFurthestPlane().normal)<epsilon)){
            allpointsonplane=0;
            break;
        }
    }
    std::cout << " All points on plane: " << allpointsonplane << std::endl;
    
 /*   std::vector<glm::vec3> verticespolytope;
    std::vector<int> indicespolytope;
    std::vector<glm::vec3> normalspolytope;
    
    
    const glm::vec3 centroidfurthestplane = getFurthestPlane().centroid;
    for (int n=0;n<farplanepoints.size()-1;n++){
        glm::vec3 pos1 = farplanepoints[n];
        glm::vec3 pos2 = farplanepoints[n+1];
        
        //add the centroid
        verticespolytope.push_back(centroidfurthestplane);
        //add as well the two vertices next two each other
        verticespolytope.push_back(pos1);
        verticespolytope.push_back(pos2);
        
        //add the normal of the triangle
        glm::vec3 norm = glm::cross(pos1-centroidfurthestplane, pos2-centroidfurthestplane);
        normalspolytope.push_back(norm);
        normalspolytope.push_back(norm);
        normalspolytope.push_back(norm);
        
    }
    
    //add the last point and the one at the begining
    
    glm::vec3 pos1 = farplanepoints[farplanepoints.size()-1];
    glm::vec3 pos2 = farplanepoints[0];
    verticespolytope.push_back(centroidfurthestplane);
    verticespolytope.push_back(pos1);
    verticespolytope.push_back(pos2);
    
    //add the normal of the triangle
    glm::vec3 norm = glm::cross(pos1-centroidfurthestplane, pos2-centroidfurthestplane);
    normalspolytope.push_back(norm);
    normalspolytope.push_back(norm);
    normalspolytope.push_back(norm);
    
    //do all indices
    for (int n=0;n<verticespolytope.size();n++){
        indicespolytope.push_back(n);
    }
    
    exportScene("fragment1.stl", verticespolytope, normalspolytope, indicespolytope);
    */
    
    
    

    
    //get the position of the centre
/*    glm::vec4 cp(getFurthestPlane().centroid.x,getFurthestPlane().centroid.y,getFurthestPlane().centroid.z,1);
    glm::mat4 trans = getTransformationForPolygoninXYPlane(getFurthestPlane());
    glm::vec4 newcentroid = cp*trans;
    farplanepointstodraw.push_back(newcentroid);
 //       std::cout <<  "centroid: " <<   newcentroid.z  << std::endl;
    
    farplanepointstodraw1.push_back(getFurthestPlane().centroid);
    
    Polygon_2    farpolygon;
    Polygon_list partition_polys;
    
    
    for (int n=0;n<farplanepoints.size();n++){
        
        glm::vec4 p(farplanepoints[n].x,farplanepoints[n].y,farplanepoints[n].z,1);
        std::cout <<  p.x << " " <<  p.y << " "  << p.z << " "   << std::endl;

      //  glm::vec4 newp = p*trans;
       // farplanepointstodraw.push_back(newp);
        //farplanepointstodraw1.push_back(farplanepoints[n]);
        //farpolygon.push_back(Point_2(newp.x, newp.y));

       // std::cout <<  newp.x << " " <<  newp.y << " "  << newp.z << " "   << std::endl;
        
    }
    
    
     std::cout <<  "&&&&&&&"   << std::endl;

    for (int n=0;n<farplanepoints.size();n++) {
        
        glm::vec4 p(farplanepoints[n].x,farplanepoints[n].y,farplanepoints[n].z,1);
      //  std::cout <<  p.x << " " <<  p.y << " "  << p.z << " "   << std::endl;
        
        glm::vec4 newp = p*trans;
        farplanepointstodraw.push_back(newp);
        farplanepointstodraw1.push_back(farplanepoints[n]);
        farpolygon.push_back(Point_2(newp.x, newp.y));
        
         std::cout <<  newp.x << " " <<  newp.y << " "  << newp.z << " "   << std::endl;
        
    }
    
    farplanepointstodraw.push_back(farplanepointstodraw[1]);
    farplanepointstodraw1.push_back(farplanepoints[1]);
    
    std::cout <<  "size: " <<   farplanepoints.size()  << std::endl;

   
    for (std::vector<glm::vec3>::reverse_iterator it = farplanepoints.rbegin();
         it != farplanepoints.rend(); ++it ) {
        glm::vec3 pos = *it;
     //   farpolygon.push_back(Point_2(pos.x, pos.y));

        
    }
    
    // check if the polygon is simple.
    std::cout << "The polygon is " <<
    (farpolygon.is_simple() ? "" : "not ") << "simple." << std::endl;
    // check if the polygon is convex
    std::cout << "The polygon is " <<
    (farpolygon.is_convex() ? "" : "not ") << "convex." << std::endl;
    
    std::cout << "The polygon is oriented " ;
    if (farpolygon.orientation() == CGAL::CLOCKWISE) std::cout <<  "Clockwise " << std::endl;
    if (farpolygon.orientation() == CGAL::COUNTERCLOCKWISE) std::cout <<  "Clockwise "<< std::endl;
    
    
    //export points
    Point_set point_set;
    
    for (VertexIterator vi = farpolygon.vertices_begin(); vi != farpolygon.vertices_end(); ++vi){
       // std::cout << "vertex " << n++ << " = " << *vi << std::endl;
       // std::cout << "vertex " << n++ << " = " << *vi->x() << std::endl;
        Point_2 &pp = *vi;
        
      //  Vertex &vert = *vi;

        
        //*vi.x();
        //Polygon_2 pp = *vi;
        point_set.insert(Point(double(pp.x()), double(pp.y()),double(0)));
    }
    

    
    /*std::cout <<  "size cgal: " <<   farpolygon.size()  << std::endl;
    if (!farpolygon.is_simple())
    {
    CGAL::approx_convex_partition_2(farpolygon.vertices_begin(),
                                    farpolygon.vertices_end(),
                                    std::back_inserter(partition_polys));
        
    }
   
    assert(CGAL::convex_partition_is_valid_2(farpolygon.vertices_begin(),
                                             farpolygon.vertices_end(),
                                             partition_polys.begin(),
                                             partition_polys.end()));*/
    
   
}

void Fragment::createSTLwithlargecones(int counterfile){
    static const double closestvalue = 0.001;
    static const double furtherstvalue = 4;
    
    std::cout << std::endl;
    std::cout << "------------pFile: " << counterfile << std::endl;
    
    //view->addGeometry(tmpqueue.front());
    //create the planes for intersection, one near the centre, the other far away from the radius of the sphere
    createPlanes(getVerticesD(), closestvalue, furtherstvalue);
    
    //std::cout << "**** " << glm::degrees(acos(glm::dot(maxdistancecentroidfragmentvertex,plane.centroid))) << std::endl;
    
    std::vector<glm::vec3> verticespolytope;
    std::vector<int> indicespolytope;
    std::vector<glm::vec3> normalspolytope;
    

    const glm::vec3 centroidfurthestplane = getFurthestPlane().centroid;
    const glm::vec3 centroidclosestplane = getClosestPlane().centroid;
    //**********************top part**********************
    for (int n=0;n<vertices.size()-1;n++){
       glm::vec3 pos1 = vertices[n];
        glm::vec3 pos2 = vertices[n+1];
        
        //add the centroid
        verticespolytope.push_back(centroidfurthestplane);
        //add as well the two vertices next two each other
        verticespolytope.push_back(pos1);
        verticespolytope.push_back(pos2);
        
        //add the normal of the triangle
        glm::vec3 norm = glm::cross(centroidfurthestplane-pos1, centroidfurthestplane-pos2);
        normalspolytope.push_back(norm);
        normalspolytope.push_back(norm);
        normalspolytope.push_back(norm);

    }
    
    //add the last point and the one at the begining
    glm::vec3 pos1 = vertices[vertices.size()-1];
    glm::vec3 pos2 = vertices[0];
    verticespolytope.push_back(centroidfurthestplane);
    verticespolytope.push_back(pos1);
    verticespolytope.push_back(pos2);
    
    //add the normal of the triangle
    glm::vec3 norm = glm::cross(pos1-centroidfurthestplane, pos2-centroidfurthestplane);
    normalspolytope.push_back(norm);
    normalspolytope.push_back(norm);
    normalspolytope.push_back(norm);

  
    
/*    verticespolytope.push_back({1,0,0});
    verticespolytope.push_back({0,0,0});
    verticespolytope.push_back({0,1,0});
    
    normalspolytope.push_back({0,0,1});
    normalspolytope.push_back({0,0,1});
    normalspolytope.push_back({0,0,1});*/
    
    //**********************bottom part**********************
   for (int n=0;n<vertices.size()-1;n++){

        glm::vec3 pos1 = vertices[n];
        glm::vec3 pos2 = vertices[n+1];
       
        //std::cout << "pos1: " << str(pos1) << std::endl;
        //std::cout << "pos2: " << str(pos2) << std::endl;

        //add the centroid
        //add as well the two vertices next two each other
        verticespolytope.push_back(pos1);
        verticespolytope.push_back({0,0,0});
        verticespolytope.push_back(pos2);
       
       
       //add the normal of the triangle
       glm::vec3 norm = glm::cross(pos2, pos1);
       
       //std::cout << "Normal: " << str(norm) << std::endl;
       normalspolytope.push_back(norm);
       normalspolytope.push_back(norm);
       normalspolytope.push_back(norm);
       
    }

    //add the last point and the one at the begining
    glm::vec3 posi1 = vertices[vertices.size()-1];
    glm::vec3 posi2 = vertices[0];
    verticespolytope.push_back(posi1);
    verticespolytope.push_back({0,0,0});
    verticespolytope.push_back(posi2);
    
    //add the normal of the triangle
    glm::vec3 normi = -glm::cross(posi2, posi1);
    normalspolytope.push_back(normi);
    normalspolytope.push_back(normi);
    normalspolytope.push_back(normi);
   



    // Writing result in STL format
    std::string filename = "..//openscad//fragmentsphere//fragment_"+std::to_string(counterfile)+".stl";

    //exportScene(filename, verticespolytope, normalspolytope, indicespolytope);
    myOwnExportSceneSTL(filename, verticespolytope, normalspolytope);
    
}
