function b = spherePolyInsideTest(poly, point)
if size(point,2) > 1
    b = zeros(1, size(point,2), 'logical');
    for i=1:numel(b)
        b(i) = spherePolyInsideTest(poly, point(:,i));
    end
else
    if 0
        % [ STILL UNDER DEVELOPMENT ]
        
        % New method, this time relying on Steorographic Project / Riemann
        % Sphere, putting the query point at the singularity. Due to the
        % conformality of the mapping, angles are preserved, and the sign
        % of the angle sum should tell me whether I am inside or outside.
        % That said, my hope is that I can avoid first applying the mapping
        % and then compute angles, as this might lead to numerical problems
        % near the singularity. Instead, I am trying to find a more direct
        % formula.
        
        point = normalize(point);  % paranoia
        poly = normalize(poly);  % paranoia
        polyNext = poly(:,[2:end,1]);
        t1 = polyNext - poly;
        t2 = normalize(cross(poly, t1));
        t1 = cross(t2, poly);
       
        if 0
            % mis-classifications if polygon vertices enclose less than pi/2:
            phiPoint = atan2(dot(repmat(point,1,size(t2,2)),t2), dot(repmat(point,1,size(t1,2)),t1));
            phiEdge = atan2(dot(polyNext,t2), dot(polyNext,t1));
        else
            if 0
                % remapping atan2 to be within [0..2pi] leads to even more
                % mis-classifications:
                phiPoint = mod(atan2(dot(repmat(point,1,size(t2,2)),t2), dot(repmat(point,1,size(t1,2)),t1)), 2*pi);
                phiEdge = mod(atan2(dot(polyNext,t2), dot(polyNext,t1)), 2*pi);
            else
                % better so far? -- stray classifications only appear once 
                % I re-enable edgeDist (else branch) again:
                phiPoint = mod(atan2(dot(repmat(point,1,size(t2,2)),t2), dot(repmat(point,1,size(t1,2)),t1)), 2*pi);
                phiEdge = mod(atan2(dot(polyNext,t2), dot(polyNext,t1)) + pi, 2*pi);
            end
        end
        positiveSide = phiPoint < phiEdge;
        
        edgeDists = zeros(1,size(poly,2));
        startDists = zeros(1,size(poly,2));
        for i=1:size(poly,2)
            [ edgeDists(i), startDists(i) ] = sphereEdgeSignedDistance(poly(:,i), polyNext(:,i), point);
        end
        dists = edgeDists;  % TODO: remove (i) once done debugging
        dists(isnan(dists)) = startDists(isnan(dists));
        [ ~, idx ] = min(abs(dists));
        idx = idx(1);
        if isnan(edgeDists(idx))
            b = positiveSide(idx);    % TODO: move t1/t2/phiXxx/positiveSide computation into this if branch
        else
            b = false;%edgeDists(idx) > 0;
        end
    else
        % Known bug: works in principle, HOWEVER, points near the polygon
        % receive (signed) winding number 0 and hence tend to be
        % misclassified as outside in 50% of the cases. Leaving it as is
        % for now, seeing that in our applications false negatives are
        % tolerable.
        %
        point = normalize(point);  % paranoia
        poly = poly - point.*sum(point.*poly, 1);   % project into plane orthogonal to point
        % poly = normalize(poly);  % not even needed, as dot and cross are bilinear forms
        c = dot(poly, poly(:, [2:end,1]));
        s = dot(repmat(point, 1, size(poly,2)), cross(poly, poly(:, [2:end,1])));
        if any(isnan(floor(0.5 + sum(atan2(s, c))./(2*pi))))
            keyboard;
        end
        b = floor(0.5 + sum(atan2(s, c))./(2*pi)) > 0;
    end
end
