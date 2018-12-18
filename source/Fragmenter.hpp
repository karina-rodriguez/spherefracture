#pragma once
#include <iostream>
#include <queue>
#include "View.hpp"
#include "Fragment.hpp"


#include <CGAL/Cartesian_d.h>
#include <cstdlib>
#include <CGAL/Min_sphere_annulus_d_traits_d.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
typedef CGAL::Cartesian_d<double>              K;
typedef CGAL::Min_sphere_annulus_d_traits_d<K> Traits;
typedef CGAL::Min_sphere_d<Traits>             Min_sphere;
typedef K::Point_d                             Point;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  EK;
typedef CGAL::Polyhedron_3<EK>                     Polyhedron_3;
typedef EK::Point_3                                Point_3;



class Fragmenter {

public:
    struct RandomFractureOptions {
        int    m;           // initial vertex count; must be m>=3
        float  niter;       // number of fractal iterations
        float  jitter;      // relative jitter of midpoint along curve
        float  amplitude;   // relative amplitude
        float  decay;       // relative amplitude/jitter decay per iteration
    };

    static RandomFractureOptions  defaultRandomFractureOptions;
    static std::vector<std::vector<glm::dvec3>> allPoints;

private:
    int radius;
    int numparts;
    int actualparts;
    View* view;

    glm::vec3 color;
    double densityline;
    double densitysphere;
    double maxpeak;


    static const double  epsilon;
    static const int evalPoints;

    static bool  epsilonSame(const glm::dvec3 &a, const glm::dvec3 &b, double epsilonScale=1.0);
    
    //set the fragment and its level
    std::queue<Fragment*> fragments;
public:

    Fragmenter(int numparts, glm::vec3 color, double radius, double densityline, double densitysphere, double maxpeak, View* view);
    Fragmenter();
    ~Fragmenter();
    int fragment();
    int testIntersections(std::vector<glm::dvec3> jline);
    int testFragment(std::vector<glm::dvec3> jline);
    void listAllFragments();
    int createPolytope();
        
    

    static bool  computeIntersection(glm::dvec3 poly1p1, glm::dvec3 poly1p2, glm::dvec3 poly2p1, glm::dvec3 poly2p2,glm::dvec3& intersectionpoint);
    static bool checkFragmentSizeSuitable(const std::vector<glm::dvec3> poly);
    static bool checkAtLeastPointsHitFragment(const int numpoints, std::vector<std::vector<glm::dvec3>> &points, const std::vector<glm::dvec3> poly);

    
    // Tim's implementation:
    bool  tryCut(const std::vector<glm::dvec3> &fragment,
                 const std::vector<glm::dvec3> &fracture,
                 std::vector<glm::dvec3> &result1,
                 std::vector<glm::dvec3> &result2);
    static std::vector<glm::dvec3>  spherePolyRandomFracture(const RandomFractureOptions &opt=defaultRandomFractureOptions);
    static bool  spherePolyIntersect(const std::vector<glm::dvec3> &poly1,
                                     const std::vector<glm::dvec3> &poly2,
                                     std::vector<glm::dvec3> &result);
    static double  spherePolyArea(const std::vector<glm::dvec3> &poly);
    static double  spherePolyAngle(const std::vector<glm::dvec3> &poly, int idx);

    static bool  spherePolyInsideTest(const std::vector<glm::dvec3> &poly, const glm::dvec3 &point);

    static bool  tests();
};
