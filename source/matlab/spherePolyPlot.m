function  spherePolyPlot(poly, colChar)
if nargin < 2
    colChar = 'k';
end
hold on;
sphereArcsPlot(poly(1,:), poly(2,:), poly(3,:), [ colChar, '--' ]);
sphereArcsPlot(poly(1,[end,1]), poly(2,[end,1]), poly(3,[end,1]), [ colChar, ':' ]);
axis equal;

if 1
    rpts = 2.*rand(3,ceil(100000/spherePolyArea(poly)))-1;
    rpts = normalize(rpts(:,dot(rpts,rpts) < 1));
    hold on;
    b = spherePolyInsideTest(poly, rpts);
    plot3(rpts(1,b), rpts(2,b), rpts(3,b), [colChar, '.']);
end

if 0
    cpoly = spherePolyConvexHull(poly);
    hold on;
    sphereArcsPlot(cpoly(1,[1:end,1]), cpoly(2,[1:end,1]), cpoly(3,[1:end,1]), [ colChar, '-' ]);
end
hold off;
