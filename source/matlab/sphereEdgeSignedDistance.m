function [ edgeDist, startDist ] = sphereEdgeSignedDistance(edgep1, edgep2, point)
if size(point,2) > 1
    edgeDist = zeros(1,size(point,2));
    startDist = zeros(1,size(point,2));
    for i=1:numel(edgeDist)
        [ edgeDist(i). startDist(i) ] = sphereEdgeSignedDistance(edgep1, edgep2, point);
    end
else
    point = normalize(point);  % paranoia
    edgep1 = normalize(edgep1);  % paranoia
    edgep2 = normalize(edgep2);  % paranoia
    A = [ edgep1, edgep2, cross(edgep1, edgep2) ];
    dt = det(A);
    edgeDist = NaN;
    startDist = NaN;
    if dt > 1e-6
        uvw = A \ point;  % u,v are used to determine quadrant spanned by p1/p2; w for side
        if all(uvw(1:2) > 0)
            % the great-circle arc that is
            % closest, so we return the signed distance to that one
            edgeDist = sign(uvw(3)) * 0.01; % * acos(dot(point, normalize(A(:,1:2) * uvw(1:2))));
        end
        startDist = acos(dot(edgep1, point));
    end
end
