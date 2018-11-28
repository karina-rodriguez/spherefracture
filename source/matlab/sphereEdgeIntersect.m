% KNOWN BUG: creates NaN if one of the spherical arcs is between directly
% opposing points. I don't see an easy fix, as there is no single smallest great-circle arc
% between the two.
% Luckily, as long as we sufficiently jitter all points, this case will not
% occur.
%
function  [ intersectionpoint, intersects ] = sphereEdgeIntersect(poly1p1, poly1p2, poly2p1, poly2p2)
epsilon = 1e-6;
intersectionpoint = [];
intersects = false;
if 0
    %
    %  Tim's implementation, using alternative derivation
    %
    edges = {};
    edges{1} = [ poly1p1, poly1p2 ];
    edges{2} = [ poly2p1, poly2p2 ];
    
    if 1
        % check for well-formed inputs
        for i=1:2
            %edges{i} = normalize(edges{i});  % shouldn't be required
            if norm(edges{i}(:,2) - edges{i}(:,1)) <= epsilon
                keyboard;
            end
        end
        if abs(dot(cross(edges{1}(:,1),edges{1}(:,2)), cross(edges{2}(:,1),edges{2}(:,2)))) >= 1 - epsilon
            warning('parallel/coplanar edges');
            keyboard;
            return;
        end
    end
    
    uvl = {};
    p = {};
    inside = [false, false];
    for i=1:2
        % Solve system of equations that traces the following path, if
        % edges are p1->p2 and q1->q2, respectively:
        %
        %   u*p1 + v*p2 - lambda*(q2-q1) - q1 = 0       <=>
        %   u*p1 + v*p2 + lambda*(q1-q2)      = q1
        %
        A = [ edges{i}, edges{3-i}(:,1)-edges{3-i}(:,2) ];
        uvl{i} = A \ edges{3-i}(:,1);
        
        inside(i) = (uvl{i}(3) > 0.5*epsilon) && (uvl{i}(3) < 1+epsilon);  % [TW:] changed to -epsilon to 0.5*epsilon
    end
    %fprintf(1, 'lambda/inside:  (%g,%g) -> (%d,%d)\n', uvl{1}(3), uvl{2}(3), inside(1), inside(2));
    %cat(2, uvl{:})
    if any(isnan(cat(1, uvl{:})))
        keyboard;
    end
    
    if all(inside)
        for i=1:2
            p{i} = uvl{i}(1).*edges{i}(:,1) + uvl{i}(2).*edges{i}(:,2);
        end
        if dot(p{1}, p{2}) > 0
            intersectionpoint = normalize(p{1});
            intersects = true;
        end
    end
else
    %
    %  Karina's implementation, using derivation jointly developed (see scribbled notes PDF)
    %
    lambda = {};
    lambda{1} = -(dot(cross(poly2p2, poly2p1), poly1p1) / dot(cross(poly2p2, poly2p1), poly1p2 - poly1p1));
    lambda{2} = -(dot(cross(poly1p2, poly1p1), poly2p1) / dot(cross(poly1p2, poly1p1), poly2p2 - poly2p1));
    
    inside = [false, false];
    for i=1:2
        inside(i) = (lambda{i} > 0.5*epsilon) && (lambda{i} < 1+epsilon);  % [TW:] changed to -epsilon to 0.5*epsilon
    end
    fprintf(1, 'lambda/inside:  (%g,%g) -> (%d,%d)\n', lambda{1}, lambda{2}, inside(1), inside(2));
    if isnan(lambda{1}) || isnan(lambda{2})
        keyboard;
    end
    intersectionpoint = [];
    if all(inside)
        intersectionpoint1 = poly1p1 + lambda{1} * (poly1p2 - poly1p1);
        intersectionpoint2 = poly2p1 + lambda{2} * (poly2p2 - poly2p1);  % [TW:] fixed multiplier to be lambda{2}, not lambda{1}
        
        %check whether these two points are actually the same
        % valid = dot(c_1, c_2) > 0
        % if the dot product is >0 then they are on the same side ==> the same
        if dot(intersectionpoint1, intersectionpoint2) > 0
            intersectionpoint = normalize(intersectionpoint1);
            intersects = true;
        end
    end
end
if 1 && intersects
    if ...
            abs(dot(intersectionpoint, cross(poly1p1, poly1p2))) > 1e-8 || ...
            abs(dot(intersectionpoint, cross(poly2p1, poly2p2))) > 1e-8
        keyboard;
    end
end
