#include "Fragmenter.hpp"

#include <assert.h>

const double  Fragmenter::epsilon = 1e-6;
const int Fragmenter::evalPoints = 2;

Fragmenter::RandomFractureOptions  Fragmenter::defaultRandomFractureOptions = {
#if 0
    3, 7, 0.28, 0.28/3.0, 0.9
#else
    4, 6, 0.4, 0.4/2.5, 0.9
#endif
};


//pointspart1= [[-100,-200,100],[-100,-200,200],[-50,-300,100],[-20,-200,150],[-200,-200,100],[-130,-200,50]];
//pointspart21= [[-300,-50,-100],[-300,-67,-300],[-300,-100,-50],[-300,-150,-250],[-900,-700,-300]];
//pointspart22 = [[-150,70,-100],[-300,0,-300],[-500,100,-250],[-300,100,-400],[-200,150,-300]];
//pointspart31=[[0,200,-50],[-120,550,-50],[-250,400,-200],[-200,400,100],[-100,500,-300]];
//pointspart32=[[0,800,200],[-200,600,0],[150,600,0],[400,700,100],[300,700,-300],[0,700,-300]];
//pointspart4=[[50,-180,50],[45,-200,-50],[-80,-300,0],[0,-300,0],[-40,-250,-50]];
//pointspart5=[[400,0,500], [400,-350,550], [400,200,500], [400,0,250], [400,150,250], [400,-200,250]];
//pointspart6=[[250,-300,-100],[250,-150,40],[120,-150,40],[190,-150,-20],[190,-150,-80]];
//pointspart7=[[-500,50,600], [-800,200,700],[-600,-180,600], [-800,0,600], [-600,-100,500]];
//pointspart8=[[-500,-50,100],[-500,200,200],[-400,300,100],[-500,100,00],[-800,400,-100],[-800,800,-100]];
//pointspart9=[[500,200,0],[250,0,0],[150,-50,-100],[150,10,-150],[250,150,-150],[300,0,-100]];
//pointspart10=[[100,300,400], [-300,400,600], [0,500,400], [-100,500,500], [200,300,300], [-400,900,700]];

static std::vector<glm::vec3> p1 ={glm::vec3(-100,-200,100),glm::vec3(-100,-200,200),glm::vec3(-50,-300,100),glm::vec3(-20,-200,150),glm::vec3(-200,-200,100),glm::vec3(-130,-200,50)};
static std::vector<glm::vec3> p2 ={glm::vec3(-300,-50,-100),glm::vec3(-300,-67,-300),glm::vec3(-300,-100,-5),glm::vec3(-300,-150,-250),glm::vec3(-900,-700,-300)};
static std::vector<glm::vec3> p3 ={glm::vec3(-150,70,-100),glm::vec3(-300,0,-300),glm::vec3(-500,100,-250),glm::vec3(-300,100,-400),glm::vec3(-200,150,-300)};
static std::vector<glm::vec3> p4 ={glm::vec3(0,200,-50),glm::vec3(-120,550,-50),glm::vec3(-250,400,-200),glm::vec3(-200,400,100),glm::vec3(-100,500,-300)};
static std::vector<glm::vec3> p5 ={glm::vec3(0,800,200),glm::vec3(-200,600,0),glm::vec3(150,600,0),glm::vec3(400,700,100),glm::vec3(300,700,-300),glm::vec3(0,700,-300)};
static std::vector<glm::vec3> p6 ={glm::vec3(50,-180,50),glm::vec3(45,-200,-50),glm::vec3(-80,-300,0),glm::vec3(0,-300,0),glm::vec3(-40,-250,-50)};
static std::vector<glm::vec3> p7 ={glm::vec3(400,0,500), glm::vec3(400,-350,550), glm::vec3(400,200,500), glm::vec3(400,0,250), glm::vec3(400,150,250), glm::vec3(400,-200,250)};
static std::vector<glm::vec3> p8 ={glm::vec3(250,-300,-100),glm::vec3(250,-150,40),glm::vec3(120,-150,40),glm::vec3(190,-150,-20),glm::vec3(190,-150,-80)};
static std::vector<glm::vec3> p9 ={glm::vec3(-500,50,600), glm::vec3(-800,200,700),glm::vec3(-600,-180,600), glm::vec3(-800,0,600), glm::vec3(-600,-100,500)};
static std::vector<glm::vec3> p10 ={glm::vec3(-500,-50,100),glm::vec3(-500,200,200),glm::vec3(-400,300,100),glm::vec3(-500,100,00),glm::vec3(-800,400,-100),glm::vec3(-800,800,-100)};
static std::vector<glm::vec3> p11 ={glm::vec3(500,200,0),glm::vec3(250,0,0),glm::vec3(150,-50,-100),glm::vec3(150,10,-150),glm::vec3(250,150,-150),glm::vec3(300,0,-100)};
static std::vector<glm::vec3> p12 ={glm::vec3(100,300,400), glm::vec3(-300,400,600), glm::vec3(0,500,400), glm::vec3(-100,500,500), glm::vec3(200,300,300), glm::vec3(-400,900,700)};

std::vector<std::vector<glm::vec3>> Fragmenter::allPoints = {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12};


/*std::vector<Fragmenter::PointsMagnets> Fragmenter::allPoints = {
    p1
   // {glm::vec3(1,0,0), glm::vec3(1,0,0)}
};*/


//////////////////////////////////////////////////////////////////////////////
// Helpers
//

template<typename T>
std::string  str(const T &val);

template<>
std::string  str(const glm::vec3 &val) {
    char  s[1024];
    sprintf(s, "(%g, %g, %g)", val.x, val.y, val.z);
    return std::string(s);
}

template<typename V>
std::string  str(const std::vector<V> &v) {
    std::string  s;
    s.append("[ ");
    for (int i=0; i<(int)v.size(); ++i) {
        if (i > 0)
            s.append(", ");
        s.append(str(v[i]));
    }
    s.append(" ]");
    return s;
}

inline double  randSigned() { return 2.0*rand()/RAND_MAX-1.0; }

inline glm::dvec3  randDirection3()
{
    glm::dvec3  dir;
    double  len;
    do {
         dir = glm::dvec3(randSigned(), randSigned(), randSigned());
         len = glm::length(dir);
    } while (len > 1);
    return (1.0/len) * dir;
}

inline glm::dmat3  randRotationMatrix3()
{
    auto  u = randDirection3();
    auto  v = glm::normalize(glm::cross(randDirection3(), u));
    auto  w = glm::cross(u, v);
    return glm::dmat3(u,v,w);
}

//////////////////////////////////////////////////////////////////////////////
// Class members
//


Fragmenter::Fragmenter(int numparts, glm::vec3 color, double radius, double densityline, double densitysphere, double maxpeak, View* view): numparts(numparts), color(color), radius(radius), densityline(densityline), densitysphere(densitysphere), maxpeak(maxpeak), view(view){
    actualparts = 0;
    //auto  jline = new JaggedLine(view->getVertexArrayID(),color,GL_LINE_LOOP,GEO_PATH, spherePolyRandomFracture());
    std::vector<glm::vec3> fracture = spherePolyRandomFracture();
    //seed two initial parts to fragment
    Fragment* fragment0 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,1.0,0.0),GL_LINE_STRIP,GEO_FRAGMENT,fracture);
    fragment0->calculateBoundingBox();
   // view->addGeometry(fragment0);
    fragments.push(fragment0);
    actualparts++;
    
    //get the vertices reversed
    std::vector<glm::vec3> thevertices = fracture;
    std::reverse(thevertices.begin(),thevertices.end());
    
    Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,1.0,0.0),GL_LINE_STRIP,GEO_FRAGMENT,thevertices);
    fragment1->calculateBoundingBox();
//    view->addGeometry(fragment1);
    fragments.push(fragment1);
    actualparts++;
    

    
   
}
Fragmenter::Fragmenter(){
    
}
Fragmenter::~Fragmenter() {
}

int Fragmenter::fragment(){

    while(actualparts<numparts){
   
        std::vector<glm::vec3> result1,result2;

        
    

        //try this cut and if it works add the two parts
        if( tryCut(fragments.front()->getVertices(), spherePolyRandomFracture(), result1,result2)){
            //remove the last fragment on the queue as we have now dealt with this
            

            fragments.pop();
            //display the fragments on screen
            glm::vec3 colorran = glm::vec3((float)rand()/RAND_MAX,(float)rand()/RAND_MAX,(float)rand()/RAND_MAX);


            //create the two new fragments
            Fragment* newfrag1 = new Fragment(view->getVertexArrayID(),colorran,GL_LINE_LOOP,GEO_FRAGMENT, result1);
            fragments.push(newfrag1);
            


            Fragment* newfrag2 = new Fragment(view->getVertexArrayID(),colorran,GL_LINE_LOOP,GEO_FRAGMENT, result2);
            fragments.push(newfrag2);
            actualparts+=1;
            

        }
        
    }

    listAllFragments();

    createPolytope();


    
    return 1;
}

int Fragmenter::testIntersections(std::vector<glm::vec3> jline){
   
    std::vector<glm::vec3> polys[2];
    polys[0] = fragments.front()->getVertices();
    polys[1] = jline;
    
    int curPoly = 0; // then one you’re “sitting on” while intersecting with the other
    int idx[2] = {0,0};
    
    int sizepoly1 = polys[0].size();
    
    //std::cout <<  jline->getNumSteps() << std::endl;
    // compute intersections
    // to be replaced for (;;) {
    for (int i=0;i<sizepoly1;i++){
        //for (;;){
        //std::cout << "ITERATION:  " << i <<  std::endl;
        
        std::cout << idx[curPoly] << std::endl;
        std::cout << (idx[curPoly] + 1) % polys[curPoly].size() << std::endl;
        
        glm::vec3 poly1p1 = polys[curPoly].at(idx[curPoly]);
        glm::vec3 poly1p2 = polys[curPoly].at((idx[curPoly] + 1) % polys[curPoly].size());
        
        
        std::cout << "Vertex actual p1: " << poly1p1.x << " " << poly1p1.y << " " << poly1p1.z << std::endl;
        std::cout << "Vertex next p1:" << poly1p2.x << " " << poly1p2.y << " " << poly1p2.z << std::endl;
        
        //iterate over the other polygon
        //in notes: for (...) { } -> idx[1-curPoly] says which edge of the other polygon intersects
        glm::vec3 intersectionpoint;
        
        for (int n=0;n<polys[1-curPoly].size();n++){
                std::cout << idx[1-curPoly] << std::endl;
                std::cout << (idx[1-curPoly] + 1) % polys[1-curPoly].size() << std::endl;
                
                glm::vec3 poly2p1 = polys[1-curPoly].at(idx[1-curPoly]);
                glm::vec3 poly2p2 = polys[1-curPoly].at((idx[1-curPoly] + 1) % polys[1-curPoly].size());
                
                std::cout << "Vertex actual p2: " << poly2p1.x << " " << poly2p1.y << " " << poly2p1.z << std::endl;
                std::cout << "Vertex next p2: " << poly2p2.x << " " << poly2p2.y << " " << poly2p2.z << std::endl;
                
                 //generate geometry to see the planes
                /* std::vector<glm::vec3> tfrag1;
                 std::vector<glm::vec3> tfrag2;
                 tfrag1.push_back(poly1p1);
                 tfrag1.push_back(poly1p2);
                 tfrag1.push_back(glm::vec3(0,0,0));
                 
                 tfrag2.push_back(poly2p1);
                 tfrag2.push_back(poly2p2);
                 tfrag2.push_back(glm::vec3(0,0,0));
                 
                 Fragment* tfragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(0.0,0.0,0.0),GL_LINE_LOOP,GEO_INTERSECTION,tfrag1);
                 Fragment* tfragment2 = new Fragment(view->getVertexArrayID(), glm::vec3(1.0,0.0,1.0),GL_LINE_LOOP,GEO_INTERSECTION,tfrag2);
                 tfragment1->calculateBoundingBox();
                 view->addGeometry(tfragment1);
                 tfragment2->calculateBoundingBox();
                 view->addGeometry(tfragment2);*/
            
                
                //calculate next point before check, so that if it intersects it does not do anything anymore ...
                idx[1-curPoly] = (idx[1-curPoly] + 1) % polys[1-curPoly].size();
                
                //COMPUTE INTERSECTION
                //if no intersection: push the current point into result and advance one
                if
                (computeIntersection(poly1p1, poly1p2, poly2p1, poly2p2,  intersectionpoint)){
                
                   // std::cout << " I INTERSECT! " << std::endl;
                    std::vector<glm::vec3> frag1;
                    std::vector<glm::vec3> frag2;
                    std::vector<glm::vec3> frag3;
                    
                    frag1.push_back(poly1p1);
                    frag1.push_back(poly1p2);
                    
                    frag2.push_back(poly2p1);
                    frag2.push_back(poly2p2);
                    
                    frag3.push_back(intersectionpoint);
                    
                    glm::vec3 colorran = glm::vec3((float)rand()/RAND_MAX,(float)rand()/RAND_MAX,(float)rand()/RAND_MAX);
                    
                    Fragment* fragment1 = new Fragment(view->getVertexArrayID(), colorran,GL_LINES,GEO_INTERSECTION,frag1);
                    fragment1->calculateBoundingBox();
                    view->addGeometry(fragment1);
                    
                    Fragment* fragment2 = new Fragment(view->getVertexArrayID(), colorran,GL_LINE_LOOP,GEO_INTERSECTION,frag2);
                    fragment2->calculateBoundingBox();
                    view->addGeometry(fragment2);
                    
                    Fragment* fragment3 = new Fragment(view->getVertexArrayID(), glm::vec3(1.0,0.0,0.0),GL_POINTS,GEO_INTERSECTION,frag3);
                    fragment3->calculateBoundingBox();
                    view->addGeometry(fragment3);
                    
                }
                
                
                
            
        }
        
        idx[curPoly] = (idx[curPoly] + 1) % polys[curPoly].size();

        
        
    }
    return 1;
}

int Fragmenter::testFragment(std::vector<glm::vec3> jline)
{
    std::vector<glm::vec3> result;
    std::vector<glm::vec3> result2;

   // std::cout << "TESTS " << std::endl;
    tryCut(fragments.front()->getVertices(), jline, result,result2);

    if (result.size()) {
        Fragment* fragment1 = new Fragment(view->getVertexArrayID(), glm::vec3(1.0,0.0,0.0),GL_LINE_STRIP,GEO_PATH,result);
        fragment1->calculateBoundingBox();
        view->addGeometry(fragment1);
        return 1;
    }
    return 0;
}

int Fragmenter::createPolytope(){
    std::queue<Fragment*> tmpqueue = fragments;
    
    int counterfile=0;
    while (!tmpqueue.empty())
    {
        Fragment* fragment = tmpqueue.front();
        fragment->createSTLwithlargecones(counterfile++);
        tmpqueue.pop();

    }
    return 1;
    

}

void Fragmenter::listAllFragments(){
    std::queue<Fragment*> tmpqueue = fragments;
    while (!tmpqueue.empty())
    {
        view->addGeometry(tmpqueue.front());

        tmpqueue.pop();
        
    }

    
}

bool  Fragmenter::computeIntersection(glm::vec3 poly1p1, glm::vec3 poly1p2, glm::vec3 poly2p1, glm::vec3 poly2p2, glm::vec3& intersectionpoint)
{
    double lambda[2];
    lambda[0] = -(glm::dot(glm::cross(poly2p2, poly2p1),poly1p1) /
                  glm::dot(glm::cross(poly2p2,poly2p1),poly1p2-poly1p1));
 
    lambda[1] = -(glm::dot(glm::cross(poly1p2,poly1p1),poly2p1) /
                  glm::dot(glm::cross(poly1p2,poly1p1),poly2p2-poly2p1));
    
    //std::cout << str(poly1p1) << " -- " << str(poly1p2) << ";   " << str(poly1p1) << " -- " << str(poly1p2) << std::endl;
    //std::cout << "LAMBDA1: "<< lambda[0] << std::endl;
    //std::cout << "LAMBDA2: "<< lambda[1] << std::endl;

    std::vector<bool> inside(2);
    for (int i=0; i<2; i++) {
        inside[i] = (lambda[i] > 0.5*epsilon) && (lambda[i] < 1+epsilon);  // [TW:] changed to -epsilon to 0.5*epsilon
    }
    if (inside[0] && inside[1]) {
        
        //compute first intersection point
        glm::vec3 vecintersectionpoint1 = glm::vec3(poly1p2-poly1p1);
        glm::vec3 vecintersectionpoint1bylambda = glm::vec3(vecintersectionpoint1.x*lambda[0], vecintersectionpoint1.y*lambda[0], vecintersectionpoint1.z*lambda[0]);
        glm::vec3 intersectionpoint1 = poly1p1 + vecintersectionpoint1bylambda;
     
        //std::cout << "Intersection point 1: "<< intersectionpoint1.x << ", " << intersectionpoint1.y << ", " << intersectionpoint1.z << std::endl;

        //compute second intersection point
        glm::vec3 vecintersectionpoint2 = glm::vec3(poly2p2-poly2p1);
        glm::vec3 vecintersectionpoint2bylambda = glm::vec3(vecintersectionpoint2.x*lambda[1], vecintersectionpoint2.y*lambda[1], vecintersectionpoint2.z*lambda[1]);
        glm::vec3 intersectionpoint2 = poly2p1 + vecintersectionpoint2bylambda;
        
        //std::cout << "Intersection point 2: "<< intersectionpoint2.x << ", " << intersectionpoint2.y << ", " << intersectionpoint2.z << std::endl;
 
        //check whether these two points are actually the same
        // valid = dot(c_1, c_2) > 0
        // if the dot product is >0 then they are on the same side ==> the same
        if (glm::dot(intersectionpoint1,intersectionpoint2)>0) {
        
            intersectionpoint = glm::normalize(intersectionpoint1);
            //std::cout << "normalised Intersection point : "<< intersectionpoint.x << ", " << intersectionpoint.y << ", " << intersectionpoint.z << std::endl;

          

            return true;
        } else
            return false;
    }

    return false;
}


bool  Fragmenter::tests()
{
    bool  success = true;
    
    std::vector<glm::vec3>  foo = { {1,0,0}, {sqrt(2.0),sqrt(2.0),0}, {0,0,1} };
    
    double                  fooArea = spherePolyArea(foo);
    //std::cout << "AREA: " << fooArea << std::endl;
    success = success && (fabs(fooArea - 0.25*M_PI) < 1e-6);

    return success;
}

bool  Fragmenter::tryCut(const std::vector<glm::vec3> &fragment,
                         const std::vector<glm::vec3> &fracture,
                         std::vector<glm::vec3> &result1,
                         std::vector<glm::vec3> &result2)
{
    //std::cout << "FRAGMENT: " << str(fragment) << std::endl;
    //std::cout << "FRACTURE: " << str(fracture) << std::endl;
    
    
    
        
    double  originalArea = spherePolyArea(fragment);
    std::cout << "originalArea " << originalArea << std::endl;

    if (spherePolyIntersect(fragment, fracture, result1)) {
        

        
        std::vector<glm::vec3>  revFracture(fracture.rbegin(), fracture.rend());
        if (spherePolyIntersect(fragment, revFracture, result2)) {
            

            
            double  area1 = spherePolyArea(result1);
            std::cout << "area1 " << area1 << std::endl;

            double  area2 = spherePolyArea(result2);
            std::cout << "area2 " << area2 << std::endl;

            double  totalArea = area1 + area2;
            std::cout << "totalArea " << totalArea << std::endl;

            double  relAreaErr = 2.0 * fabs(totalArea - originalArea) / (totalArea + originalArea);
            std::cout << "relAreaErr " << relAreaErr << std::endl;

                if (relAreaErr < 1e-5 &&                                       // single pieces, please
                             std::max(area1 / area2, area2 / area1) < 4.0)  // maximum area ratio
            
                    if
                    ((fragments.size()>1) && (checkFragmentSizeSuitable(result1)&&checkFragmentSizeSuitable(result2))) // pieces of radius less than 0.9
          //  std::cout << " ---->> " << checkAtLeastPointsHitFragment(1,allPoints,result1) << " --- " << checkAtLeastPointsHitFragment(1,allPoints,result2) << std::endl;
                    return
                        /*(evalPoints) &&*/
                        (checkAtLeastPointsHitFragment(1,allPoints,result1)||checkAtLeastPointsHitFragment(1,allPoints,result2));
            
        }
    }
    return false;
    
}

// spherePolyRandomFracture produces spherical polygon of m * 2^niter points.
std::vector<glm::vec3>  Fragmenter::spherePolyRandomFracture(const RandomFractureOptions &opt)
{
    auto  A = randRotationMatrix3();

    std::vector<glm::vec3>  poly;
    {
        poly.reserve(opt.m);
        const double  f = 2.0 * M_PI / opt.m;
        for (int i=0; i<opt.m; ++i) {
            glm::vec3  v((float)(f*i + f*opt.jitter * randSigned()),
                         (float)(f*opt.amplitude * randSigned()),
                         0.f);
            poly.push_back(A * glm::vec3(cos(v.x)*cos(v.y), sin(v.x)*cos(v.y), sin(v.y)));
        }
    }

    {
        double  jitter = opt.jitter;
        double  amplitude = opt.amplitude;
        for (int iter=0; iter<opt.niter; ++iter) {
            std::vector<glm::vec3>  newPoly;
            newPoly.reserve(2 * poly.size());
            for (int idx=0; idx<(int)poly.size(); ++idx) {
                int  idxSucc = (idx + 1) % ((int)poly.size());

                auto  w = glm::normalize(glm::cross(poly[idx], poly[idxSucc]));
                auto  u = poly[idx];
                auto  v = glm::cross(w, u);
    
                double  arclen = acos(glm::dot(poly[idxSucc], u));
                double  midarclen = arclen * (0.5+jitter * randSigned());
                double  midampl = arclen * opt.amplitude * randSigned();
                auto  midpoint =
                    u * float(cos(midarclen) * cos(midampl))
                    + v * float(sin(midarclen) * cos(midampl))
                    + w * float(sin(midampl));

                newPoly.push_back(poly[idx]);
                newPoly.push_back(midpoint);
            }
            poly.swap(newPoly);
        
            jitter *= opt.decay;
            amplitude *= opt.decay;
        }
    }
    
    return poly;
}

bool  Fragmenter::epsilonSame(const glm::vec3 &a, const glm::vec3 &b, double epsilonScale)
{
    return glm::length(a - b) <= epsilonScale*epsilon;
}

//! Intersection of two spherical polygons.
/*! Polygons' orientation is clockwise around their inner region.
 *
 *  Returns true if input polygons actually intersect AND if the resulting polygon is non-empty.
 *
 *  If intersection would yield multiple components, \a result
 *  contains only one of them, and computeIntersection() still returns
 *  true. (The case of multiple components is explicitly tested for by
 *  in tryCut(), which relies on that behaviour.)
 * 
 *  NOTE: this code heavily relies on computeIntersection to *exclude*
 *  half(!) an epsilon ball of each line around its starting point and
 *  to *include* a full epsilon ball around its end point. I already
 *  adjusted computeIntersection accordingly.
 */
bool  Fragmenter::spherePolyIntersect(const std::vector<glm::vec3> &poly1,
                                      const std::vector<glm::vec3> &poly2,
                                      std::vector<glm::vec3> &result)
{
    const std::vector<glm::vec3>  *polys[2] = {&poly1, &poly2};
#if 0
    // for debug purposes, reverse cut polygon:
    std::vector<glm::vec3>  revPoly2(poly2.rbegin(), poly2.rend());
    polys[1] = &revPoly2;
#endif

    if (!tests()) {
        std::cerr << "TESTS FAILED.\n";
        exit(1);
    }
    
    int  curPoly = 0; // then one you’re “sitting on” while intersecting with the other
    int  idx = 0;

    result.clear();
    result.push_back((*polys[curPoly])[idx]);
    // idx now points to the "previous vertex" on "current polygon";
    // sometimes, this is not the same as the last entry in result,
    // when an intersection point has been added.
    int  idxSucc = (idx + 1) % (*polys[curPoly]).size();  // index of end point of current edge on current polygon
    
    bool  noCut = true;
    
    for (int steps=0; steps<1000; steps++) {
        bool       intersects = false;
        int        jdx=-1, jdxSucc=-1;
        glm::vec3  intersection;
        for (jdx=0; jdx<(int)(*polys[1-curPoly]).size(); jdx++) {
            jdxSucc = (jdx+1) % (*polys[1-curPoly]).size();
            if (computeIntersection(result.back(), (*polys[curPoly])[idxSucc],
                                    (*polys[1-curPoly])[jdx], (*polys[1-curPoly])[jdxSucc],
                                    intersection)) {
                intersects = true;
                break;
            }
        }
        //std::cout << "intersects: " << intersects << " " << jdx << "/" << (int)(*polys[1-curPoly]).size() << std::endl;
        
        if (!intersects) {
            // no intersection, so we add end point of this edge and
            // proceed on current polygon:
            idx = idxSucc;
            idxSucc = (idx + 1) % (*polys[curPoly]).size();
            if (!epsilonSame(result.back(), (*polys[curPoly])[idx], 2.0))
                result.push_back((*polys[curPoly])[idx]);
        } else {
            noCut = false;
            if (glm::dot(result.back(),
                         glm::cross((*polys[1-curPoly])[jdx], (*polys[1-curPoly])[jdxSucc])) < 0) {
                // we are approaching from outside the cut polygon ->
                // discard what has been collected so far and start
                // from intersection point:
                result.clear();
                result.push_back(intersection);
                // leave idx as is, even if "behind" intersection point
            } else {
                // we are leaving the cut polygon, so add intersection
                // point and switch polygon:
                assert(result.size());
                if (!epsilonSame(result.back(), intersection, 2.0))
                    result.push_back(intersection);
                idx = jdx;
                idxSucc = jdxSucc;
                curPoly = 1 - curPoly;
                // leave idx as is, even if "behind" intersection point
            }
        }

        if (result.size() > 1 && epsilonSame(result.front(), result.back(), 2.0)) {
            result.pop_back();
            break;
        }

        // [TW:] it looks like this termination criterion was poorly motivated and even introduces a bug:
        //if (curPoly == 0 && idx == 0)
        //    break;
    }

    //std::cout << str(result) << std::endl;
    return !noCut && result.size() > 0;
}

double  Fragmenter::spherePolyArea(const std::vector<glm::vec3> &poly)
{
    double  sum = 0;
    for (int i=0; i<(int)poly.size(); ++i)
        sum += spherePolyAngle(poly, i);
    return sum - M_PI * ((int)poly.size() - 2);
}

double  Fragmenter::spherePolyAngle(const std::vector<glm::vec3> &poly, int idx)
{
    int  n = (int)poly.size();
    glm::vec3  pp = poly[(idx + n -1) % n];   // previous
    glm::vec3  pc = poly[idx];                // current
    glm::vec3  pn = poly[(idx + 1) % n];      // next
#if 1  // Paranoia
    pp = glm::normalize(pp);
    pc = glm::normalize(pc);
    pn = glm::normalize(pn);
#endif
    glm::vec3  np = glm::normalize(glm::cross(pp, pc));
    glm::vec3  nn = glm::normalize(glm::cross(pc, pn));

    glm::vec3  tp = glm::cross(np, pc);

    return fmod(atan2(glm::dot(nn, tp), glm::dot(nn, np)) + M_PI, 2.0*M_PI);
}

// UNTESTED port from Matlab...
// KNOWN BUG (as observed in Matlab): works in principle, HOWEVER,
// for particularly jaggy polygons, points near the polygon receive
// (signed) winding number 0 and hence tend to be misclassified as
// outside in 50% of the cases. Leaving it as is for now, seeing that
// in our applications false negatives are tolerable.
//
bool  Fragmenter::spherePolyInsideTest(const std::vector<glm::vec3> &poly, const glm::vec3 &point)
{
    auto  pl = poly;                   // will be modified

    auto  pt = glm::normalize(point);  // paranoia
    
    for (int i=0; i<(int)pl.size(); ++i)
        pl[i] = pl[i] - dot(pt, pl[i]) * pt;  // project into plane orthogonal to pt
    
    double  sum = 0.0;
    for (int i=0; i<(int)pl.size(); ++i) {
        int  iSucc = (i+1) % (int)pl.size();
        sum += atan2(glm::dot(pt, glm::cross(pl[i], pl[iSucc])), glm::dot(pl[i], pl[iSucc]));
    }
    
   std::cout << "###### " << floor(0.5 + sum / (2*M_PI)) << std::endl;
    return floor(0.5 + sum / (2*M_PI)) > 0;
}

bool Fragmenter::checkFragmentSizeSuitable(const std::vector<glm::vec3> poly){
    
    
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^"  << std::endl;

    //min sphere for
    Point  P[poly.size()];
    double coord[3];
    
    // Build also the points for the convex hull
    std::vector<Point_3> points;

    for (int n=0;n<poly.size();n++){
      //  std::cout << "Point: " << (double)poly[n].x << ", " << (double)poly[n].y << ", " << (double)poly[n].z << std::endl;
        
        coord[0] = (double)poly[n].x;
        coord[1] = (double)poly[n].y;
        coord[2] = (double)poly[n].z;
        P[n] = Point(3, coord,coord+3);
        
        int dimension = CGAL::Feature_dimension<Point, K>::value;
        assert(dimension == 0);
        
        points.push_back(Point_3(poly[n].x,poly[n].y,poly[n].z));

        
    }
    
    Min_sphere  ms (P, P+3);             // smallest enclosing sphere
    double  radius = sqrt(ms.squared_radius());
    std::cout << "   radius: " << radius << " " ;
    //   CGAL::set_pretty_mode (std::cout);
    //   std::cout << ms1;
    
    
    //calculate the convex hull
    Polyhedron_3 polyhedron;
    
    // compute convex hull of non-collinear points
    CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);
    
    
    CGAL::Side_of_triangle_mesh<Polyhedron_3, EK> inside(polyhedron);
  //  int nb_inside = 0;
  //  CGAL::Bounded_side res = inside(Point_3(0,0,0));
  //  if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    std::cout << (inside(Point_3(0,0,0))==CGAL::ON_BOUNDED_SIDE) << std::endl;
    
    
    return
    sqrt(ms.squared_radius())<0.9 &&      //check the fragment is not too big
        (inside(Point_3(0,0,0))!=CGAL::ON_BOUNDED_SIDE);     //check the origin is not inside the fragment
    
}
bool Fragmenter::checkAtLeastPointsHitFragment(const int numpoints,const std::vector<std::vector<glm::vec3>> points,const std::vector<glm::vec3> poly){

    for (int i=0;i<points.size();i++){
        
        std::vector<glm::vec3> pointsforFragment = points[i];
        int numfound=0;
        for (int j=0;j<pointsforFragment.size();j++){
            std::cout << "**pos to check:***>>> " << pointsforFragment[j].x << ",  " << pointsforFragment[j].y << ", " << pointsforFragment[j].z ;

            if (spherePolyInsideTest(poly, pointsforFragment[j])) numfound++;
                if (numfound==numpoints) {
                    std::cout << "******>>> " << i << " and " << j << std::endl;
                    numfound++;
                    removePointsForThisFragment(numfound);
                    return 1;
                }
            
        }
        
    }
   
    return 0;
}
                
void Fragmenter::removePointsForThisFragment(const int elem){
    
    std::remove(allPoints.begin(), allPoints.end(), allPoints[elem]);
    
    
}

