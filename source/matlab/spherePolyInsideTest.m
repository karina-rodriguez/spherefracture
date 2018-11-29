% Known bug: works in principle, HOWEVER, points near the polygon
% receive (signed) winding number 0 and hence tend to be
% misclassified as outside in 50% of the cases. Leaving it as is 
% for now, seeing that in our applications false negatives are
% tolerable.
%
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
    c = dot(poly, poly(:, [2:end,1]));
    s = dot(repmat(point, 1, size(poly,2)), cross(poly, poly(:, [2:end,1])));
    b = floor(0.5 + sum(atan2(s, c))./(2*pi)) > 0;
end
