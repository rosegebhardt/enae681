% Flow parameters
alpha = 0; % angle of attack (rad)
Gamma_v = 1; % vorticity (rad/sec)
delta0 =  -0.2; % thickness parameter with zero camber (unitless)
beta = 0.0; % camber parameter (unitless)
r0 = 1; % cylinder radius (m)
u0 = 1; % freestream velocity (m/s)
k0 = r0*(delta0 + 1); % uncambered thickness (unitless)
rho = 1000; % water density (kg/m^3)

% Derived parameters
A0 = pi*r0^2*(1-((delta0+1)^4)/((1-delta0^2)^2));
deltaConstant = (sqrt(pi*r0^2-A0)-sqrt(pi*r0^2))/...
                (sqrt(pi*r0^2-A0)+sqrt(pi*r0^2));

% Circle and foil shapes definition
delta = deltaConstant*sqrt(1-beta^2);
x0 = r0*(delta + 1i*beta);
k = r0*(delta + sqrt(1 - beta^2));

% Export parameters as a .mat file
save parameterData.mat alpha Gamma_v delta0 r0 u0 k0 rho beta delta x0 k