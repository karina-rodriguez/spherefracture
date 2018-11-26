#include "Fragmenter.hpp"
#include "Fragment.hpp"
#include "JaggedLine.hpp"

Fragmenter::Fragmenter(int numparts, glm::vec3 color, float radius, float densityline, float densitysphere, float maxpeak, View* view): numparts(numparts), color(color), radius(radius), densityline(densityline), densitysphere(densitysphere), maxpeak(maxpeak), view(view){
    actualparts = 0;
    actuallevel = 0;
    JaggedLine* jline = new JaggedLine(view->getVertexArrayID(),color,GL_LINE_LOOP,GEO_PATH, densityline,maxpeak, true,glm::vec3(0,1,0));
    
    //seed two initial parts to fragment
    Fragment* fragment0 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,1.0,0.0),GL_LINE_STRIP,GEO_FRAGMENT,jline->getVertices());
    fragment0->calculateBoundingBox();
    fragment0->setLevel(0);
    view->addGeometry(fragment0);
    fragments.push_back(fragment0);
    actualparts++;
    
    Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3((float)rand()/RAND_MAX,(float)rand()/RAND_MAX,(float)rand()/RAND_MAX),GL_LINE_LOOP,GEO_FRAGMENT,jline->getVerticesReversed());
    fragment1->calculateBoundingBox();
    fragment1->setLevel(0);
//    view->addGeometry(fragment1);
    fragments.push_back(fragment1);
    actualparts++;
    
}
Fragmenter::~Fragmenter() {
}


int Fragmenter::fragment(){
   


    
    
    //while (actualparts<numparts){
    
    
    //Go through all our fragments
    //commented for the moment
/*    bool foundFragment =false;
    while (!fragments.empty())
    {
        Fragment* fragment = fragments.back();
        //get a fragment from our list which has not been processed
        if (fragment->getLevel()==actuallevel){
            foundFragment = true;
            //remove it from our list as we are going to process it
            fragments.pop_back();
            if(!generateTwoFragments(fragment,jline)) {
                std::cout << "something went wrong cutting" << std::endl;
                break;
            }
        }
        if (!foundFragment) actuallevel++;
    }
    
    listAllFragments();*/
    

    return 1;
    
}
    
bool Fragmenter::generateTwoFragments(Fragment* fragment, JaggedLine* jline){

    // for (;;) {
   //     for (...) { } -> idx[1-curPoly] says which edge of the other polygon intersects
    //}

    /*std::vector<Vector3> polys[2]; // each is a spherical polygon
    int curPoly = 0; // then one you’re “sitting on” while intersecting with the other
    int idx[2] = {0,0};
    std::vector<Vector3> result;
    for (;;) {
        // compute intersections
        for (...) { } -> idx[1-curPoly] says which edge of the other polygon intersects
        if no intersection:
            idx[curPoly] = (idx[curPoly] + 1) % poly[curPoly].size();
        else:
            if intersection takes wrong turn:
                result.clear();
        result.push_back(intersection point);
            else:
                if new intersection point is different enough from poly[1-curPoly]:
                    result.push_back(intersection point);
                else
                    result.push_back(poly[1-curPoly][idx[1-curPoly]]);
        curPoly = 1 - curPoly;
        idx[curPoly] = (idx[curPoly] + 1) % poly[curPoly].size();
        if result.size() >= 2 && length(result.front() - result.back()) < epsilon)
            break;
    }
    result.pop_back();*/
    

}

int Fragmenter::testFragment(JaggedLine* jline){
    
    std::vector<glm::vec3> polys[2];
    polys[0] = fragments.front()->getVertices();
    polys[1] = jline->getVertices();
    
    int curPoly = 0; // then one you’re “sitting on” while intersecting with the other
    int idx[2] = {0,0};
    std::vector<glm::vec3> result;
    
    
    int index =0;
    int sizepoly1 = polys[0].size();
    std::vector<glm::vec3>::iterator itr0 = polys[0].begin();
    glm::vec3 polyinitial = *itr0;
    
    //std::cout <<  jline->getNumSteps() << std::endl;
    // compute intersections
    // to be replaced for (;;) {
    for (int i=0;i<jline->getNumSteps();i++){
    //for (;;){
        //std::cout << "ITERATION:  " << i <<  std::endl;

        //std::cout << idx[curPoly] << std::endl;
        //std::cout << (idx[curPoly] + 1) % polys[curPoly].size() << std::endl;
        
        glm::vec3 poly1p1 = polys[curPoly].at(idx[curPoly]);
        glm::vec3 poly1p2 = polys[curPoly].at((idx[curPoly] + 1) % polys[curPoly].size());
        
        //std::cout << "Vertex actual p1: " << poly1p1.x << " " << poly1p1.y << " " << poly1p1.z << std::endl;
        //std::cout << "Vertex next p1:" << poly1p2.x << " " << poly1p2.y << " " << poly1p2.z << std::endl;
        
        //iterate over the other polygon
        //in notes: for (...) { } -> idx[1-curPoly] says which edge of the other polygon intersects
        bool intersect = 0;
        glm::vec3 intersectionpoint;

        for (int n=0;n<polys[1-curPoly].size();n++){
            if (!intersect){
                //std::cout << idx[1-curPoly] << std::endl;
                //std::cout << (idx[1-curPoly] + 1) % polys[1-curPoly].size() << std::endl;
                
                glm::vec3 poly2p1 = polys[1-curPoly].at(idx[1-curPoly]);
                glm::vec3 poly2p2 = polys[1-curPoly].at((idx[1-curPoly] + 1) % polys[1-curPoly].size());
                
                //std::cout << "Vertex actual p2: " << poly2p1.x << " " << poly2p1.y << " " << poly2p1.z << std::endl;
                //std::cout << "Vertex next p2: " << poly2p2.x << " " << poly2p2.y << " " << poly2p2.z << std::endl;
                
                /* generate geometry to see the planes
                 std::vector<glm::vec3> frag1;
                std::vector<glm::vec3> frag2;
                frag1.push_back(poly1p1);
                frag1.push_back(poly1p2);
                frag1.push_back(glm::vec3(0,0,0));
                
                frag2.push_back(poly2p1);
                frag2.push_back(poly2p2);
                frag2.push_back(glm::vec3(0,0,0));
                
                Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,0.0,0.0),GL_LINE_LOOP,GEO_INTERSECTION,frag1);
                Fragment* fragment2 = new Fragment(view->getVertexArrayID(), glm::vec3(1.0,0.0,1.0),GL_LINE_LOOP,GEO_INTERSECTION,frag2);
                fragment1->calculateBoundingBox();
                view->addGeometry(fragment1);
                fragment2->calculateBoundingBox();
                view->addGeometry(fragment2);*/

                
                //calculate next point before check, so that if it intersects it does not do anything anymore ...
                idx[1-curPoly] = (idx[1-curPoly] + 1) % polys[1-curPoly].size();
                
                //COMPUTE INTERSECTION
                //if no intersection: push the current point into result and advance one
                if (computeIntersection(poly1p1,poly1p2,poly2p1,poly2p2,intersectionpoint)&&(!intersect)){
                    
                    //TO DO....
                    // checkRightTurn(poly1p1,poly1p2,poly2p1,poly2p2);
                    intersect=1;
                    
                }

                
            
            }
        }
        // if nothing intersect, then push back the point
        if (!intersect) {
            //std::cout << "no intersection, then insert same point "<< idx[curPoly] << std::endl;
            //std::cout << "Vertex inserted: " << poly1p1.x << " " << poly1p1.y << " " << poly1p1.z << std::endl;

            result.push_back(polys[curPoly].at(idx[curPoly]));
            
        }
        //if intersect then calculate the point
        else{
            //insert initial point of this vector
            result.push_back(polys[curPoly].at(idx[curPoly]));

            curPoly = 1 - curPoly;
            //std::cout << "normalised Intersection point : "<< intersectionpoint.x << ", " << intersectionpoint.y << ", " << intersectionpoint.z << std::endl;
            result.push_back(intersectionpoint);




        }
        //and increase our index for the current polygon whoever it is
        
       // std::cout << "num "<< (idx[curPoly] + 1)  << std::endl;
       // std::cout << "num "<< (idx[curPoly] + 1)  << std::endl;

        idx[curPoly] = (idx[curPoly] + 1) % polys[curPoly].size();
        if ((result.size() >= 2) && (length(result.front() - result.back()) < epsilon)) {
            break;
        }
        
        
        
    }
    
    Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(1.0,0.0,0.0),GL_LINE_STRIP,GEO_PATH,result);
    fragment1->calculateBoundingBox();
    view->addGeometry(fragment1);
}
void Fragmenter::listAllFragments(){
    std::cout << "**********Fragments: " << fragments.size() << std::endl;
    for (std::vector<Fragment*>::iterator it = fragments.begin() ; it != fragments.end(); ++it){
        Fragment* frag = *it;
        std::cout << "Vertices " << frag->getVertices().size() << " " << frag->getLevel() << std::endl;
        
    }

}
bool Fragmenter::computeIntersection(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2, glm::vec3& intersectionpoint){
    double lambda[2];
    lambda[0] = -glm::dot(glm::cross(poly2p2, poly2p1),poly1p1)/
    glm::dot(glm::cross(poly2p2,poly2p1),poly1p2-poly1p1);
 
    lambda[1] = -glm::dot(glm::cross(poly1p2,poly1p1),poly2p1)/
    glm::dot(glm::cross(poly1p2,poly1p1),poly2p2-poly2p1);
    
    
    //std::cout << "LAMBDA1: "<< lambda[0] << std::endl;
    //std::cout << "LAMBDA2: "<< lambda[1] << std::endl;

    std::vector<bool> inside(2);
    for (int i=0; i<2; i++){
        inside[i] = (lambda[i] > -epsilon) && (lambda[i] < 1.000006);
    }
    if (inside[0]&&inside[1]) {
        
        //compute first intersection point
        glm::vec3 vecintersectionpoint1 = glm::vec3(poly1p2-poly1p1);
        glm::vec3 vecintersectionpoint1bylambda = glm::vec3(vecintersectionpoint1.x*lambda[0], vecintersectionpoint1.y*lambda[0], vecintersectionpoint1.z*lambda[0]);
        glm::vec3 intersectionpoint1 = poly1p1 + vecintersectionpoint1bylambda;
     
        //std::cout << "Intersection point 1: "<< intersectionpoint1.x << ", " << intersectionpoint1.y << ", " << intersectionpoint1.z << std::endl;

        //compute second intersection point
        glm::vec3 vecintersectionpoint2 = glm::vec3(poly2p2-poly2p1);
        glm::vec3 vecintersectionpoint2bylambda = glm::vec3(vecintersectionpoint2.x*lambda[0], vecintersectionpoint2.y*lambda[0], vecintersectionpoint2.z*lambda[0]);
        glm::vec3 intersectionpoint2 = poly2p1 + vecintersectionpoint2bylambda;
        
        //std::cout << "Intersection point 2: "<< intersectionpoint2.x << ", " << intersectionpoint2.y << ", " << intersectionpoint2.z << std::endl;
 
        //check whether these two points are actually the same
        // valid = dot(c_1, c_2) > 0
        // if the dot product is 0 then they are the same
        if (glm::dot(intersectionpoint1,intersectionpoint2)>0){
        
            intersectionpoint = glm::normalize(intersectionpoint1);
           // std::cout << "normalised Intersection point : "<< intersectionpoint.x << ", " << intersectionpoint.y << ", " << intersectionpoint.z << std::endl;

            //create two geometries to see the lines
            /*std::vector<glm::vec3> frag1;
            std::vector<glm::vec3> frag2;
            frag1.push_back(poly1p1);
            frag1.push_back(poly1p2);
            frag1.push_back(glm::vec3(0,0,0));
            
            frag2.push_back(poly2p1);
            frag2.push_back(poly2p2);
            frag2.push_back(glm::vec3(0,0,0));
            
            Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,1.0,0.0),GL_LINE_LOOP,GEO_INTERSECTION,frag1);
            Fragment* fragment2 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,0.0,1.0),GL_LINE_LOOP,GEO_INTERSECTION,frag2);
            fragment1->calculateBoundingBox();
            view->addGeometry(fragment1);
            fragment2->calculateBoundingBox();
            view->addGeometry(fragment2);*/

            return true;
        }else return false;
    }

    return false;
}
bool Fragmenter::checkRightTurn(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2){
    
    
}
