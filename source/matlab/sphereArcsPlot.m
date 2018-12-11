function  sphereArcsPlot(vert_x, vert_y, vert_z, style)

if 0
    plot3(vert_x, vert_y, vert_z, style);
else
    step = pi/180;
    xyz = zeros(3,0);
    for i=1:numel(vert_x)-1
        v1 = [ vert_x(i); vert_y(i); vert_z(i) ];
        v2 = [ vert_x(i+1); vert_y(i+1); vert_z(i+1) ];
        d = norm(v2-v1);
        da = 2*asin(0.5*d);
        n = ceil(da / step);
        t = da .* ([0:n-1] ./ n - 0.5);
        t = sin(t) ./ d + 0.5;
        uvw = (1-t).*v1 + t.*v2;
        uvw = normalize(uvw);
        xyz = cat(2, xyz, uvw);
    end
    xyz(:,end+1) = [ vert_x(end); vert_y(end); vert_z(end) ];
    plot3(xyz(1,:), xyz(2,:), xyz(3,:), style);
    plot3(vert_x, vert_y, vert_z, '*k');
end
