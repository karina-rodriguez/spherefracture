#pragma once
#include <iostream>
#include "View.hpp"
#include "Fragment.hpp"

#include "JaggedLine.hpp"

class Fragmenter {
    
private:
    int radius;
    int numparts;
    int actualparts;
    int actuallevel;
    View* view;

    glm::vec3 color;
    double densityline;
    double densitysphere;
    double maxpeak;
    
    static const double  epsilon;
    static bool  epsilonSame(const glm::vec3 &a, const glm::vec3 &b, double epsilonScale=1.0);
    
    //set the fragment and its level
    std::vector<Fragment*> fragments;
public:

    Fragmenter(int numparts, glm::vec3 color, double radius, double densityline, double densitysphere, double maxpeak, View* view);
    Fragmenter();
    ~Fragmenter();
    int fragment();
    int testIntersections(JaggedLine* jline);
    int testFragment(JaggedLine* jline);
    void listAllFragments();
        
    
    static bool  computeIntersection(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2,glm::vec3& intersectionpoint);
    static bool  checkRightTurn(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2);

    // Tim's implementation:
    bool  tryCut(const std::vector<glm::vec3> &fragment,
                 const std::vector<glm::vec3> &fracture,
                 std::vector<glm::vec3> &result1,
                 std::vector<glm::vec3> &result2);
    static bool  spherePolyIntersect(const std::vector<glm::vec3> &poly1,
                                     const std::vector<glm::vec3> &poly2,
                                     std::vector<glm::vec3> &result);
    static double  spherePolyArea(const std::vector<glm::vec3> &poly);
    static double  spherePolyAngle(const std::vector<glm::vec3> &poly, int idx);

    static bool  tests();
};
