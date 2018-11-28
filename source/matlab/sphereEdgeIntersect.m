function  [ intersectionpoint, intersects ] = sphereEdgeIntersect(poly1p1, poly1p2, poly2p1, poly2p2)
epsilon = 1e-6;
lambda = {};
lambda{1} = -(dot(cross(poly2p2, poly2p1), poly1p1) / ...
    dot(cross(poly2p2, poly2p1), poly1p2 - poly1p1));
lambda{2} = -(dot(cross(poly1p2, poly1p1), poly2p1) / ...
    dot(cross(poly1p2, poly1p1), poly2p2 - poly2p1));

inside = [false, false];
for i=1:2
    inside(i) = (lambda{i} > 0.5*epsilon) && (lambda{i} < 1+epsilon);  % [TW:] changed to -epsilon to 0.5*epsilon
end

intersects = false;
intersectionpoint = [];
if all(inside)
    % compute first intersection point
    vecintersectionpoint1 = poly1p2 - poly1p1;
    vecintersectionpoint1bylambda = lambda{1} * vecintersectionpoint1;
    intersectionpoint1 = poly1p1 + vecintersectionpoint1bylambda;
    
    % compute second intersection point
    vecintersectionpoint2 = poly2p2 - poly2p1;
    vecintersectionpoint2bylambda = lambda{1} * vecintersectionpoint2;
    intersectionpoint2 = poly2p1 + vecintersectionpoint2bylambda;
    
    %check whether these two points are actually the same
    % valid = dot(c_1, c_2) > 0
    % if the dot product is >0 then they are on the same side ==> the same
    if dot(intersectionpoint1,intersectionpoint2) > 0
        intersectionpoint = normalize(intersectionpoint1);
        
        intersects = true;
    end
end
