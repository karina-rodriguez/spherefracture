#include "JaggedLine.hpp"

JaggedLine::JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, double density,double maxpeak,bool anticlockwise): Geometry(vertexarrayIDT,colour, primitive, type), density(density), maxpeak(maxpeak) {
//    const double PI = 3.141592653589793238462643383279502884197;
    //generate a random rotation by giving a random vector
    const double radius = 1;
    double randomx = ((float)rand()/RAND_MAX * 2) - 1;
    double randomy = ((float)rand()/RAND_MAX * 2) - 1;
    double randomz = ((float)rand()/RAND_MAX * 2) - 1;
    std::cout << "random : " << randomx << ", " << randomy <<", " << randomz << std::endl;
    
    glm::mat4 rotation = glm::rotate((float)((double)rand()/RAND_MAX *360),glm::vec3(randomx,randomy,randomz));
  //  glm::mat4 rotation = glm::rotate((float)90,glm::vec3(1,0,0));

    //phi is angle in XZ plane, theta is angle in XY plane
    // iterate through all aximuth angles, but only for one altitude: PI/2
    double  phi_step =  glm::pi<float>()/density;
    //check if we need to do one phi_step less to avoid having two points similar to each other
    for (double phi = 0; phi < glm::two_pi<double>()-phi_step; phi+=phi_step){// azimuth angle
        double theta = glm::pi<double>()/2;//altitude angle
        
        double randomvalue = ((double)rand()/RAND_MAX * (2*(0.6*phi_step))) - (0.6*phi_step);
        double ph = phi + randomvalue;
        
        /*
        double x;
        double xtmp = radius * cos(ph) * sin(theta);
        if (anticlockwise) x = -xtmp;
        else*/
        double x = radius * cos(ph) * sin(theta);
        double y = radius * sin(ph) * sin(theta);
        double ztmp = radius * cos(theta);
        //get a random value between -maxpeak*2 and maxpeak*2
        double randomzvalue = ((double)rand()/RAND_MAX * (maxpeak*2)) - maxpeak;
        double z = ztmp + randomzvalue;
       //  double z = radius * cos(theta);
//std::cout << "here : " << x << ", " << y <<", " << z << std::endl;
        //-x makes it go anticlockwise
        glm::vec4 p(-x,y,z,0);
   //     std::cout << "p : " << p.x << ", " << p.y <<", " << p.z << std::endl;
        glm::vec4 newp = p * rotation;
       std::cout << "newp : " << newp.x << ", " << newp.y <<", " << newp.z << std::endl;
        vertices.push_back(glm::vec3(newp.x,newp.y,newp.z));
//  vertices.push_back(glm::vec3(p.x,p.y,p.z));
    }
    
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
    
  //  std::cout << "**********$$$$$$$$$$$$$$$$" << vertices.size() << std::endl;

    
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), &colors[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
}
JaggedLine::~JaggedLine() {
}
//testing jagged line
JaggedLine::JaggedLine(GLuint vertexarrayIDT, glm::vec3 colour, GLenum primitive, geo_type type, double density,double maxpeak,bool anticlockwise, glm::vec3 rot): Geometry(vertexarrayIDT,colour, primitive, type), density(density), maxpeak(maxpeak) {
    //    const double PI = 3.141592653589793238462643383279502884197;
    //generate a random rotation by giving a random vector
    const double radius = 1;
    double randomx = ((float)rand()/RAND_MAX * 2) - 1;
    double randomy = ((float)rand()/RAND_MAX * 2) - 1;
    double randomz = ((float)rand()/RAND_MAX * 2) - 1;
    std::cout << "random : " << randomx << ", " << randomy <<", " << randomz << std::endl;
    glm::mat4 rotation = glm::rotate((float)90.0,rot);
    //  glm::mat4 rotation = glm::rotate((float)90,glm::vec3(1.0,0.0,0.0));
    
    //phi is angle in XZ plane, theta is angle in XY plane
    // iterate through all aximuth angles, but only for one altitude: PI/2
    double  phi_step =  glm::pi<double>()/density;
    //check if we need to do one phi_step less to avoid having two points similar to each other
    for (double phi = 0; phi < glm::two_pi<double>()-phi_step; phi+=phi_step){// azimuth angle
        double theta = glm::pi<double>()/2;//altitude angle
        
        double randomvalue = ((double)rand()/RAND_MAX * (2*(0.6*phi_step))) - (0.6*phi_step);
        double ph = phi;// + randomvalue;
        
        /*
         double x;
         double xtmp = radius * cos(ph) * sin(theta);
         if (anticlockwise) x = -xtmp;
         else*/
        double x = radius * cos(ph) * sin(theta);
        double y = radius * sin(ph) * sin(theta);
        double ztmp = radius * cos(theta);
        //get a random value between -maxpeak*2 and maxpeak*2
        double randomzvalue = ((double)rand()/RAND_MAX * (maxpeak*2)) - maxpeak;
        double z = ztmp;// + randomzvalue;
        //  double z = radius * cos(theta);
        //std::cout << "here : " << x << ", " << y <<", " << z << std::endl;
        //-x makes it go anticlockwise
        glm::vec4 p(-x,y,z,0);
        //     std::cout << "p : " << p.x << ", " << p.y <<", " << p.z << std::endl;
        glm::vec4 newp = p * rotation;
        std::cout << "newp : " << newp.x << ", " << newp.y <<", " << newp.z << std::endl;
        vertices.push_back(glm::vec3(newp.x,newp.y,newp.z));
        //  vertices.push_back(glm::vec3(p.x,p.y,p.z));
       // if (num==5) break;
     //   num++;
    }
    
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
    
   // std::cout << "**********$$$$$$$$$$$$$$$$" << vertices.size() << std::endl;
    
    
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec3), &colors[0], GL_STATIC_DRAW);
    
    glGenBuffers(1, &normalbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);
    
}/*
void
function fractaltest(useExistingFigure)

% Debug infrastructure
%
if nargin < 1
close all;
m = 3;
n = 4;
figure;
for i=1:m*n
subplot(m, n, i);
fractaltest(true);
end
return;
else
if ~useExistingFigure
figure;
end
end

% Parameters
%
%   two reasonable set of parameters follow
%
if 0
jitter = 0.28;              % relative jitter of midpoint along curve
amplitude = jitter / 3.0;   % relative amplitude
decay = 0.9;                % relative amplitude/jitter decay per iteration
niter = 5;                  % number of fractal iterations
m = 3;                      % initial vertex count; must be m>=3
else
jitter = 0.4;               % relative jitter of midpoint along curve
amplitude = jitter / 2.5;   % relative amplitude
decay = 0.9;                % relative amplitude/jitter decay per iteration
niter = 4;                  % number of fractal iterations
m = 4;                      % initial vertex count; must be m>=3
end

% Algorithm (produces m * 2^niter points)
%
x = (2*pi)*(0:m-1)/m;
x = x + ((jitter/m)*2*pi).*(2*rand(size(x))-1);
y = ((amplitude/m)*2*pi).*(2*rand(size(x))-1);

A = rand(3);
A(:,2) = cross(A(:,3), A(:,1));
A(:,3) = cross(A(:,1), A(:,2));
A = normalize(A);  % column-wise normalisation
%A'*A    % test

xyz = A * cat(1, cos(x).*cos(y), sin(x).*cos(y), sin(y));

phidense = linspace(0,2*pi,1000);
circ1 = A * cat(1, cos(phidense), sin(phidense), zeros(size(phidense)));
circ2 = A * cat(1, zeros(size(phidense)), cos(phidense), sin(phidense));
circ3 = A * cat(1, sin(phidense), zeros(size(phidense)), cos(phidense));

for iter=1:niter
xyz2 = xyz(:,[2:end,1]);
d = xyz2 - xyz;
len = sqrt(sum(d.^2,1));

w = normalize(cross(xyz, xyz2));
u = xyz;
v = cross(w, u);

arclen = acos(dot(xyz2, u));
midarclen = arclen .* (0.5+jitter.*(2*rand(size(arclen))-1));
midampl = arclen .* amplitude .* (2*rand(size(arclen))-1);
midpoint = u .* (cos(midarclen).*cos(midampl)) + v .* (sin(midarclen).*cos(midampl)) + w .* sin(midampl);

xyz = reshape(cat(1, xyz, midpoint), 3, []);  % interleave xyz and midpoint vector arrays

jitter = decay * jitter;
amplitude = decay * amplitude;
end

%plot([x,2*pi],y([1:end,1]),'-k');
plot3(xyz(1,:), xyz(2,:), xyz(3,:), 'b-');
axis equal;
hold on;
plot3(xyz(1,[end,1]), xyz(2,[end,1]), xyz(3,[end,1]), 'b--');
plot3(circ1(1,:), circ1(2,:), circ1(3,:), 'r:');
plot3(circ2(1,:), circ2(2,:), circ2(3,:), 'k:');
plot3(circ3(1,:), circ3(2,:), circ3(3,:), 'k:');
hold off;

if 0
% test spherePolyInsideTest()
%
%   Result so far: works in principle, HOWEVER, points near the polygon
%   receive (signed) winding number 0 and hence tend to be
%   misclassified as outside in 50% of the cases.
%
rpts = 2.*rand(3,1000)-1;
rpts = normalize(rpts(:,dot(rpts,rpts) < 1));
hold on;
b = spherePolyInsideTest(xyz, rpts);
plot3(rpts(1,b), rpts(2,b), rpts(3,b), '.g');
plot3(rpts(1,~b), rpts(2,~b), rpts(3,~b), '.r');
hold off;
end

function b = spherePolyInsideTest(poly, point)
if size(point,2) > 1
b = zeros(1, size(point,2), 'logical');
for i=1:numel(b)
b(i) = spherePolyInsideTest(poly, point(:,i));
end
else
point = normalize(point);  % paranoia
poly = poly - point.*sum(point.*poly, 1);   % project into plane orthogonal to point
poly = normalize(poly);
c = dot(poly(:, [1:end,1]), poly(:, [2:end,1,2]));
s = dot(repmat(point, 1, 1+size(poly,2)), cross(poly(:, [1:end,1]), poly(:, [2:end,1,2])));
b = floor(0.5 + sum(atan2(s, c))./(2*pi)) > 0;
end


function x = normalize(x)
x = x ./ sqrt(sum(x.^2,1))
*/
