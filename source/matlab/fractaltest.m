function fractaltest(useExistingFigure)

% Debug infrastructure
%
if nargin < 1
    close all;
    m = 2;
    n = 2;
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

% Parameters
%
if 1
    jitter = 0.28;             % relative jitter of midpoint along curve
    amplitude = jitter / 3.0;   % relative amplitude
    decay = 0.9;                % relative amplitude/jitter decay per iteration
    niter = 5;                  % number of fractal iterations
    m = 3;                      % initial vertex count; must be m>=3
else
    jitter = 0.333;             % relative jitter of midpoint along curve
    amplitude = jitter / 3.0;   % relative amplitude
    decay = 0.9;                % relative amplitude/jitter decay per iteration
    niter = 5;                  % number of fractal iterations
    m = 4;                      % initial vertex count; must be m>=3
end

% Algorithm (produces m * 2^niter points)
%
x = (2*pi)*(0:m-1)/m;
x = x + ((jitter/m)*2*pi).*(2*rand(size(x))-1);
y = ((amplitude/m)*2*pi).*(2*rand(size(x))-1);

A = rand(3);
A(:,2) = cross(A(:,3), A(:,1));
A(:,3) = cross(A(:,1), A(:,2));
A = normalize(A);  % column-wise normalisation
%A'*A    % test

xyz = A * cat(1, cos(x).*cos(y), sin(x).*cos(y), sin(y));

phidense = linspace(0,2*pi,1000);
circ1 = A * cat(1, cos(phidense), sin(phidense), zeros(size(phidense)));
circ2 = A * cat(1, zeros(size(phidense)), cos(phidense), sin(phidense));
circ3 = A * cat(1, sin(phidense), zeros(size(phidense)), cos(phidense));

for iter=1:niter
    xyz2 = xyz(:,[2:end,1]);
    d = xyz2 - xyz;
    len = sqrt(sum(d.^2,1));
    
    w = normalize(cross(xyz, xyz2));
    u = xyz;
    v = cross(w, u);
    
    arclen = acos(dot(xyz2, u));
    midarclen = arclen .* (0.5+jitter.*(2*rand(size(arclen))-1));
    midampl = arclen .* amplitude .* (2*rand(size(arclen))-1);
    midpoint = u .* (cos(midarclen).*cos(midampl)) + v .* (sin(midarclen).*cos(midampl)) + w .* sin(midampl);
    
    xyz = reshape(cat(1, xyz, midpoint), 3, []);  % interleave xyz and midpoint vector arrays
    
    jitter = decay * jitter;
    amplitude = decay * amplitude;
end

%plot([x,2*pi],y([1:end,1]),'-k');
plot3(xyz(1,:), xyz(2,:), xyz(3,:), 'b-');
axis equal;
hold on;
plot3(xyz(1,[end,1]), xyz(2,[end,1]), xyz(3,[end,1]), 'b--');
plot3(circ1(1,:), circ1(2,:), circ1(3,:), 'r:');
plot3(circ2(1,:), circ2(2,:), circ2(3,:), 'k:');
plot3(circ3(1,:), circ3(2,:), circ3(3,:), 'k:');
hold off;

function x = normalize(x)
x = x ./ sqrt(sum(x.^2,1));
