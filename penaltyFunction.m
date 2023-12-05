function penalty = penaltyFunction(zVect)

load parameterData.mat r0 x0 k

% Map vortex position to cylinder plane
zv1 = zVect(1); zv2 = zVect(2);
zv = zv1 + 1i*zv2;
xv = 0.5*(zv + sign(zv1).*(zv.^2 - 4*k^2).^0.5);

% Distance from center of cylinder (related to distance from foil)
% dist < 1 ==> foil interior, dist > 1 ==> foil exterior
distance = abs(xv-x0)/r0;
distance(distance < 1) = nan; % remove interior points

% Interior penalty function equation (penalty -> Inf on boundary)
% penalty = 2.5e-3./((distance-1).^2); % alternate penalty
penalty = 1e-2./(distance-1);
