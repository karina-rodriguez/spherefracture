function [ result1, result2, success ] = tryCut(fragment, fracture)

originalArea = spherePolyArea(fragment);
fprintf(1, 'originalArea: %g\n', originalArea);

[ result1, ~, success ] = spherePolyIntersect(fragment, fracture);
if success
    [ result2, ~, success ] = spherePolyIntersect(fragment, fracture(:,end:-1:1));
    if success
        area1 = spherePolyArea(result1);
        area2 = spherePolyArea(result2);
        totalArea = area1 + area2;
        
        relAreaErr = 2.0 * abs(totalArea - originalArea) ./ (totalArea + originalArea);
        
        success = relAreaErr < 0.001 ...              % maximum relative error of split (to catch non-simple polygons)
            && max(area1 / area2, area2 / area1) < 4  % maximum area ratio
        fprintf(1, 'a1, a2, sum, relErr; success:  %g, %g, %g, %g; %d\n', area1, area2, totalArea, relAreaErr, success);
    end
end
