function tests

%fractaltest;

pool = { spherePolyRandomFracture() };
pool{end+1} = fliplr(pool{1});
plotPool(pool);

while numel(pool) < 10
    fprintf(1, 'fragment count: %d\n', numel(pool));
    idx = ceil(numel(pool)*rand(1));
    fracture = spherePolyRandomFracture();
    [ result1, result2, success ] = tryCut(pool{idx}, fracture);
    if success
        pool{idx} = result1;
        pool{end+1} = result2;
        plotPool(pool);
    end
end

function plotPool(pool)
colChars = 'rgbcmyk';

figure;
if 1
    phidense = linspace(0,2*pi,1000);
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
