function  poly = spherePolyConvexHull(poly)
[~,i] = min(poly(1,:));  % unless polygon lies entirely in constant-x plane, this vertex will lie on the convex hull
% Okay, if the point above was /definitely/ on the hull, I could simply
% compute the hull in n steps; just in order to not risk this degenerate
% case, I'd take the risk of doubling the compute time by requiring one
% complete round with no change...
nochange = 0;
while nochange < size(poly,2)
    if spherePolyAngle(poly, i) >= pi
        poly = poly(:,[1:(i-1),(i+1):end]);
        i = 1 + mod(i-2, size(poly,2));
        nochange = 0;
    else
        nochange = nochange + 1;
        i = 1 + mod(i, size(poly,2));
    end
end
