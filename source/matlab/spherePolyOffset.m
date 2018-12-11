% [ STILL UNDER DEVELOPMENT ]
function poly = spherePolyOffset(poly, offset)
% Strategies to avoid the most obvious problematic cases:
% NO:
% - Corners with inner angle less than 90 degree get chamfered.
% - Both during chamfering and moving of vertices, no vertex 
%   shall be moved further than its projection can travel on
%   the adjacent edges without ...
% YES:
% - move only one edge at a time, and check for shift collapsing adjacent
%   edges.

polyNext = poly(:,[2:end,1]);
for i=1:size(poly,2)
    % extract local polyline of this edge and its two adjacent edges (four
    % vertices); in the case of triangles, first and last point will be
    % identical:
    pl = poly(:,1+mod((i-1:i+2)-1, size(poly,2)));
    
    if 0
        sidesBefore = onPositiveSide(pl(:,2:3), pl(:,[1,4]));
        %...
        sidesAfter = onPositiveSide(pl(:,2:3), pl(:,[1,4]));
    else
        % here's the ad-hoc and dangerous way of offsetting:
        % AND YET: there's even a bug beyond that. Since I incorporated it
        % in this one-edge-at-the-time loop, some vertices are moved
        % multiple times. I should reimplement the naive method, this time
        % moving each vertex based on the two adjacent edges:
        out2 = outward2(pl(:,2:3));
        pl(:,2:3) = normalize(pl(:,2:3) + tan(offset).*out2);
    end
    
    % TODO: once "proper" branch is finnished, it will signal point
    % deletions via setting parts of pl to NaN. This will have to be tested
    % for here:
    poly(:,1+mod((i-1:i+2)-1, size(poly,2))) = pl;
end

% returns two(!) outward directions, one for each vertex, as those
% directions differ when on the sphere
function out = outward2(edge)
dir = edge(:,2) - edge(:,1);
out = cat(2, cross(dir, edge(:,1)), cross(dir, edge(:,2)));

function b = onPositiveSide(edge, points)
dir = edge(:,2) - edge(:,1);
outward = cat(1, dir(2,:), -dir(1,:));
points = points - edge(:,1);
b = dot(points, outward) > 0;
