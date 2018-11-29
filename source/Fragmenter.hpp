#pragma once
#include <iostream>
#include "View.hpp"
#include "Fragment.hpp"

#include "JaggedLine.hpp"

class Fragmenter {

public:
    struct RandomFractureOptions {
        int    m;           // initial vertex count; must be m>=3
        float  jitter;      // relative jitter of midpoint along curve
        float  amplitude;   // relative amplitude
        float  decay;       // relative amplitude/jitter decay per iteration
        float  niter;       // number of fractal iterations
    };

    static RandomFractureOptions  defaultRandomFractureOptions;
    
private:
    int radius;
    int numparts;
    int actualparts;
    int actuallevel;
    
    glm::vec3 color;
    float densityline;
    float densitysphere;
    float maxpeak;
    View* view;
    
    static const double  epsilon;
    static bool  epsilonSame(const glm::vec3 &a, const glm::vec3 &b, double epsilonScale=1.0);
    
    //set the fragment and its level
    std::vector<Fragment*> fragments;
public:
    Fragmenter(int numparts, glm::vec3 color, float radius, float densityline, float densitysphere, float maxpeak, View* view);
    ~Fragmenter();
    int fragment();
    int testFragment(JaggedLine* jline);
    bool generateTwoFragments(Fragment* fragment, JaggedLine* jline);
    void listAllFragments();
    
    static bool  computeIntersection(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2,glm::vec3& intersectionpoint);
    static bool  checkRightTurn(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2);

    // Tim's implementation:
    bool  tryCut(const std::vector<glm::vec3> &fragment,
                 const std::vector<glm::vec3> &fracture,
                 std::vector<glm::vec3> &result1,
                 std::vector<glm::vec3> &result2);
    static std::vector<glm::vec3>  spherePolyRandomFracture(const RandomFractureOptions &opt=defaultRandomFractureOptions);
    static bool  spherePolyIntersect(const std::vector<glm::vec3> &poly1,
                                     const std::vector<glm::vec3> &poly2,
                                     std::vector<glm::vec3> &result);
    static double  spherePolyArea(const std::vector<glm::vec3> &poly);
    static double  spherePolyAngle(const std::vector<glm::vec3> &poly, int idx);

    static bool  tests();
};
