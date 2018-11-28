function area = spherePolyArea(poly)
sum = 0;
for i=1:size(poly,2)
    sum = sum + spherePolyAngle(poly, i);
end
area = sum - pi * (size(poly,2) - 2);
