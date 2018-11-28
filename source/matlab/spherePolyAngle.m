function a = spherePolyAngle(poly, idx)

n = size(poly,2);
pp = poly(:,1+mod(idx-2,n));  % previous
pc = poly(:,idx);             % current
pn = poly(:,1+mod(idx,n));     % next
if 1  % Paranoia
    pp = normalize(pp);
    pc = normalize(pc);
    pn = normalize(pn);
end
np = normalize(cross(pp, pc));
nn = normalize(cross(pc, pn));

tp = cross(np, pc);

a = mod(atan2(dot(nn, tp), dot(nn, np)) + pi, 2*pi);
