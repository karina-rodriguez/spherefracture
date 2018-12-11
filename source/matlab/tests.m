function tests

close all;
%fractaltest;

if 0
    if 1
        pool = { spherePolyRandomFracture() };
    else
        %pool = { [ 1, 0, -1, 0; 0, 1, 0, -1; 0, 0, 0, 0 ] };
        pool{1} = pool{1} + 0.01*(2*rand(size(pool{1}))-1);
    end
else
    start_poly = load('debug_start_poly_1.mat');
    pool = { start_poly.poly };
    clear start_poly;
end
if 1
    poly = pool{1};
    save('latest_start_poly.mat', 'poly');
    clear poly;
end

%pool{end+1} = fliplr(pool{1});
plotPool(pool);
return;
if 0
    %fracture = spherePolyRandomFracture();
    fracture = pool{1}([2,3,1],:);
    
    if 1
        [Q,R] = qr(rand(3));
        fracture = Q * fracture;
    end
    
    hold on;
    spherePolyPlot(fracture);
    
    pool2 = {};
    for i=1:2
        [ result1, ~, success, intersectionPoints ] = spherePolyIntersect(pool{i}, fracture);
        if success
            hold on;
            %plot3(intersectionPoints(1,:),intersectionPoints(2,:),intersectionPoints(3,:),'*m');
            %keyboard;
            pool2{end+1} = result1;
            [ result2, ~, success, intersectionPoints ] = spherePolyIntersect(pool{i}, fracture(:,end:-1:1));
            if success
                pool2{end+1} = result2;
            end
        end
    end
    
    plotPool(pool2);
    keyboard;
    return;
end

while numel(pool) < 10
    fprintf(1, 'fragment count: %d\n', numel(pool));
    idx = ceil(numel(pool)*rand(1));
    fracture = spherePolyRandomFracture();
    [ result1, result2, success ] = tryCut(pool{idx}, fracture);
    if success
        pool{idx} = result1;
        pool{end+1} = result2;
        %plotPool(pool);
    end
end

if 0
    for i=1:numel(pool)
        pool{i} = spherePolyOffset(pool{i}, -1.0 * (pi/180));
    end
end
plotPool(pool);

function plotPool(pool)
colChars = 'rgbcmyk';

figure;
if 1
    phidense = linspace(0,2*pi,10000);
    A = eye(3);
    circ1 = A * cat(1, cos(phidense), sin(phidense), zeros(size(phidense)));
    circ2 = A * cat(1, zeros(size(phidense)), cos(phidense), sin(phidense));
    circ3 = A * cat(1, sin(phidense), zeros(size(phidense)), cos(phidense));
    
    plot3(circ1(1,:), circ1(2,:), circ1(3,:), 'k:');
    axis equal;
    hold on;
    plot3(circ2(1,:), circ2(2,:), circ2(3,:), 'k:');
    plot3(circ3(1,:), circ3(2,:), circ3(3,:), 'k:');
    hold off;
end

for i=1:numel(pool)
    spherePolyPlot(pool{i}, colChars(1+mod(i-1,numel(colChars))));
end
