function poly = spherePolyRandomFracture(opt, debug_level)

if nargin < 1 || isempty(opt)
    % Parameters
    %
    %   two reasonable set of parameters follow
    %
    opt = struct();
    if 0
        opt.m = 3;                      % initial vertex count; must be m>=3
        opt.jitter = 0.28;              % relative jitter of midpoint along curve
        opt.amplitude = opt.jitter / 3.0;   % relative amplitude
        opt.decay = 0.9;                % relative amplitude/jitter decay per iteration
        opt.niter = 5;                  % number of fractal iterations
    else
        opt.m = 4;                      % initial vertex count; must be m>=3
        opt.jitter = 0.4;               % relative jitter of midpoint along curve
        opt.amplitude = opt.jitter / 2.5;   % relative amplitude
        opt.decay = 0.9;                % relative amplitude/jitter decay per iteration
        opt.niter = 4;                  % number of fractal iterations
    end
end
if nargin < 2
    debug_level = 0;
end

% Algorithm (produces m * 2^niter points)
%
x = (2*pi)*(0:opt.m-1)/opt.m;
x = x + (opt.jitter*(2*pi/opt.m)).*(2*rand(size(x))-1);
y = (opt.amplitude*(2*pi/opt.m)).*(2*rand(size(x))-1);

A = rand(3);
A(:,2) = cross(A(:,3), A(:,1));
A(:,3) = cross(A(:,1), A(:,2));
A = normalize(A);  % column-wise normalisation
%A'*A    % test

poly = A * cat(1, cos(x).*cos(y), sin(x).*cos(y), sin(y));

for iter=1:opt.niter
    polyNext = poly(:,[2:end,1]);
    d = polyNext - poly;
    
    w = normalize(cross(poly, polyNext));
    u = poly;
    v = cross(w, u);
    
    arclen = acos(dot(polyNext, u));
    midarclen = arclen .* (0.5+opt.jitter.*(2*rand(size(arclen))-1));
    midampl = arclen .* opt.amplitude .* (2*rand(size(arclen))-1);
    midpoint = u .* (cos(midarclen).*cos(midampl)) + v .* (sin(midarclen).*cos(midampl)) + w .* sin(midampl);
    
    poly = reshape(cat(1, poly, midpoint), 3, []);  % interleave poly and midpoint vector arrays
    
    opt.jitter = opt.decay * opt.jitter;
    opt.amplitude = opt.decay * opt.amplitude;
end

if debug_level > 0
    phidense = linspace(0,2*pi,1000);
    circ1 = A * cat(1, cos(phidense), sin(phidense), zeros(size(phidense)));
    circ2 = A * cat(1, zeros(size(phidense)), cos(phidense), sin(phidense));
    circ3 = A * cat(1, sin(phidense), zeros(size(phidense)), cos(phidense));
    
    %plot([x,2*pi],y([1:end,1]),'-k');
    plot3(poly(1,:), poly(2,:), poly(3,:), 'b-');
    axis equal;
    hold on;
    plot3(poly(1,[end,1]), poly(2,[end,1]), poly(3,[end,1]), 'b--');
    plot3(circ1(1,:), circ1(2,:), circ1(3,:), 'r:');
    plot3(circ2(1,:), circ2(2,:), circ2(3,:), 'k:');
    plot3(circ3(1,:), circ3(2,:), circ3(3,:), 'k:');
    hold off;
end
