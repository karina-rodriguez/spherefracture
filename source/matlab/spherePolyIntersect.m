% Intersection of two spherical polygons.
% Polygons' orientation is clockwise around their inner region.
%
%  Returns true if input polygons actually intersect AND if the resulting polygon is non-empty.
%
%  If intersection would yield multiple components, \a result
%  contains only one of them, and computeIntersection() still returns
%  true. (The case of multiple components is explicitly tested for by
%  in tryCut(), which relies on that behaviour.)
% 
%  NOTE: this code heavily relies on computeIntersection to *exclude*
%  half(!) an epsilon ball of each line around its starting point and
%  to *include* a full epsilon ball around its end point. I already
%  adjusted computeIntersection accordingly.
%
function [ result, intersected, success ] = spherePolyIntersect(poly1, poly2)

polys = {poly1, poly2};
if 0
    % for debug purposes, reverse cut polygon:
    polys{2} = fliplr(polys{2});
end

curPoly = 1; % then one you?re ?sitting on? while intersecting with the other
idx = 1;

result = zeros(3,0);
result(:,end+1) = polys{curPoly}(:,idx);
% idx now points to the "previous vertex" on "current polygon";
% sometimes, this is not the same as the last entry in result,
% when an intersection point has been added.
idxSucc = 1+mod(idx, size(polys{curPoly},2));  % index of end point of current edge on current polygon

intersected = false;

for steps=1:1000
    intersects = false;
    jdx=0;
    jdxSucc=0;
    for jdx=1:size(polys{3-curPoly},2)
        jdxSucc = 1+mod(jdx+1,size(polys{3-curPoly},2));
        [ intersection, success ] = sphereEdgeIntersect( ...
            result(:,end), polys{curPoly}(:,idxSucc), ...
            polys{3-curPoly}(:,jdx), polys{3-curPoly}(:,jdxSucc));
        if success
            intersects = true;
            break;
        end
    end
    
    if ~intersects
        % no intersection, so we add end point of this edge and
        % proceed on current polygon:
        idx = idxSucc;
        idxSucc = 1+mod(idx, size(polys{curPoly},2));
        if ~epsilonSame(result(:,end), polys{curPoly}(:,idx), 2.0);
            result(:,end+1) = polys{curPoly}(:,idx);
        end
    else
        intersected = true;
        if dot(result(:,end), cross(polys{3-curPoly}(:,jdx), polys{3-curPoly}(:,jdxSucc))) > 0
            % we are approaching from outside the cut polygon ->
            % discard what has been collected so far and start
            % from intersection point:
            result = zeros(3,0);
            result(:,end+1) = intersection;
            % leave idx as is, even if "behind" intersection point
        else
            % we are leaving the cut polygon, so add intersection
            % point and switch polygon:
            if isempty(result)
                error('never reached');
            end
            if ~epsilonSame(result(:,end), intersection, 2.0)
                result(:,end+1) = intersection;
            end
            idx = jdx;
            idxSucc = jdxSucc;
            curPoly = 3 - curPoly;
            % leave idx as is, even if "behind" intersection point
        end
    end
    
    if size(result,2) > 1 && epsilonSame(result(:,1), result(:,end), 2.0)
        result = result(:,end-1);
        break;
    end
    
    if curPoly == 1 && idx == 1
        break;
    end
end

success = intersected && ~isempty(result);


function b = epsilonSame(u, v, multiplier)
if nargin < 3
    multiplier = 1;
end
epsilon = 1e-6;
u = u - v;
b = sqrt(dot(u,u)) <= multiplier * epsilon;
