function fractaltest(useExistingFigure)

% Debug infrastructure
%
if nargin < 1
    close all;
    m = 3;
    n = 4;
    figure;
    for i=1:m*n
        subplot(m, n, i);
        fractaltest(true);
    end
    return;
else
    if ~useExistingFigure
        figure;
    end
end

poly = spherePolyRandomFracture([], 1);

if 0
    % test spherePolyInsideTest()
    %
    rpts = 2.*rand(3,1000)-1;
    rpts = normalize(rpts(:,dot(rpts,rpts) < 1));
    hold on;
    b = spherePolyInsideTest(poly, rpts);
    plot3(rpts(1,b), rpts(2,b), rpts(3,b), '.g');
    plot3(rpts(1,~b), rpts(2,~b), rpts(3,~b), '.r');
    hold off;
end

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
    c = dot(poly(:, [1:end,1]), poly(:, [2:end,1,2]));
    s = dot(repmat(point, 1, 1+size(poly,2)), cross(poly(:, [1:end,1]), poly(:, [2:end,1,2])));
    b = floor(0.5 + sum(atan2(s, c))./(2*pi)) > 0;
end
