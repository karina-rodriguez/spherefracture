#pragma once
// Include GLM
#include <glm/glm.hpp>
#include "Geometry.hpp"
#include <iostream>
#include <list>
#include <fstream>
#include "common/objloader.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPK;
typedef CGAL::Partition_traits_2<EPK>                       EPKTraits;
typedef EPKTraits::Point_2                                  Point_2;
typedef EPKTraits::Polygon_2                                Polygon_2;
typedef Polygon_2::Vertex_iterator                          Vertex_iterator;
typedef std::list<Polygon_2>                                Polygon_list;
typedef CGAL::Creator_uniform_2<int, Point_2>               Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator>   Point_generator;
typedef Polygon_2::Vertex_iterator VertexIterator;

/*


#include <CGAL/Point_set_3.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
//#include <CGAL/Point_set_3/IO.h>
#include <fstream>


//#include <limits>


//For pointset in CGAL

 typedef K::FT FT;
 typedef K::Point_3 Point;
 typedef K::Vector_3 Vector;
 typedef CGAL::cpp11::array<unsigned char, 3> Color;
 typedef CGAL::Point_set_3<Point> Point_set;
 typedef Point_set::Property_map<Color> Color_map;
 typedef Point_set::Property_map<FT> FT_map;*/

struct Plane{
    glm::dvec3 normal;
    glm::dvec3 centroid;
    
};
class Fragment : public Geometry {


    

private:
    Plane plane;
    Plane closestplane;
    Plane furthestplane;
    //stores the points of the planar polygons close and far from the centre
    std::vector<glm::dvec3> closeplanepoints;
    std::vector<glm::dvec3> farplanepoints;
    double maxdistancecentroidfragment;
    glm::dvec3 maxdistancecentroidfragmentvertex;


    static int seedOrientation;
    static const double  epsilon;

public:

    std::vector<glm::dvec3> farplanepointstodraw1;
    std::vector<glm::dvec3> farplanepointstodraw;
    Fragment(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, std::vector<glm::dvec3> verticest);
    ~Fragment();
    int createPlanes(const std::vector<glm::dvec3> vertices, const double close,  const double far);
    int checkIntersectionwithPlane(glm::dvec3 point, Plane plane, glm::dvec3& result);
    Plane getPlane();
    Plane getClosestPlane();
    Plane getFurthestPlane();
    glm::mat4 getTransformationForPolygoninXYPlane(Plane theplane);
    void createSTL(int counterfile);
    void createSTLwithlargecones(int counterfile);
    static void setSeedOrientation(const int&  seed);
    

};
