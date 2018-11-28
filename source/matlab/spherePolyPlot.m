function  spherePolyPlot(poly, colChar)
if nargin < 2
    colChar = 'k';
end
hold on;
plot3(poly(1,:), poly(2,:), poly(3,:), [ colChar, '--' ]);
plot3(poly(1,[end,1]), poly(2,[end,1]), poly(3,[end,1]), [ colChar, ':' ]);
axis equal;

rpts = 2.*rand(3,10000)-1;
rpts = normalize(rpts(:,dot(rpts,rpts) < 1));
hold on;
b = spherePolyInsideTest(poly, rpts);
plot3(rpts(1,b), rpts(2,b), rpts(3,b), [colChar, '.']);
hold off;
