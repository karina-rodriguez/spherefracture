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
    
    glm::vec3 color;
    float densityline;
    float densitysphere;
    float maxpeak;
    View* view;
    const float epsilon = 1e-6;
    
    //set the fragment and its level
    std::vector<Fragment*> fragments;
public:
    Fragmenter(int numparts, glm::vec3 color, float radius, float densityline, float densitysphere, float maxpeak, View* view);
    ~Fragmenter();
    int fragment();
    int testFragment(JaggedLine* jline);
    bool generateTwoFragments(Fragment* fragment, JaggedLine* jline);
    void listAllFragments();
    bool computeIntersection(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2,glm::vec3& intersectionpoint);
    bool checkRightTurn(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2);
        
};
