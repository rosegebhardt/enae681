function drag = dragFunction(zVect) % dragFunction(zv1,zv2)

load parameterData.mat alpha Gamma_v r0 u0 rho beta delta x0 k

% Shape definition of cylinder and foil
offsetFactor = 1.05; % small offset from foil surface, prevents integration over singularities
angle = linspace(0,2*pi,2000);
circle = offsetFactor*r0*exp(1i*angle) + x0;
foil = circle + (offsetFactor*k)^2./circle;

% Vortex location cylinder and foil domains
zv1 = zVect(1); zv2 = zVect(2);
zv = zv1 + 1i*zv2;
xv = 0.5*(zv + sign(zv1).*(zv.^2 - 4*k^2).^0.5);

% Location of image vortex (models vortex interaction with foil)
xim = x0 + r0^2./conj(xv-x0);

% Bound vortex strength needed to match Kutta condition
Gamma_0 = Gamma_v*((abs(xv-x0)^2 - r0^2)./(abs(xv-k)^2)) + ...
          4*pi*u0*((delta-k)*sin(alpha) - beta*cos(alpha));

% Surface velocity in the cylinder domain      
unmappedVelocity = u0*exp(1i*alpha) - u0*exp(-1i*alpha)*r0^2./conj((circle-x0).^2) + ...
                   1i*Gamma_0./(2*pi*conj(circle-x0)) + ...
                   1i*Gamma_v./(2*pi*conj(circle-xv)) - ...
                   1i*Gamma_v./(2*pi*conj(circle-xim));

% Surface velocity in the foil domain                 
mappedVelocity = unmappedVelocity.*conj((circle.^2)./(circle.^2  - k^2));

% Surface pressure (only finite values)
surfacePressure = -rho*abs(mappedVelocity).^2/2000; % pressure in kPa
finiteIndices = isfinite(surfacePressure);
surfacePressure = surfacePressure(finiteIndices);

% Integrate surface pressure to get forces
drag = -cumtrapz(imag(foil(finiteIndices)),surfacePressure); drag = drag(end); % force in kN
% lift = cumtrapz(real(foil(finiteIndices)),surfacePressure); lift = lift(end);  % force in kN
