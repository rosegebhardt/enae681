%%  Generate visuals for optimization task

% SECTION 1: generates sample images of the flow field and the pressure
% field to illustrate the physical system.

% SECTION 2: calculates the value of the drag function and penalty function
% over the complex domain.

% SECTION 3: produces a contour plot of the drag and penalty functions over
% the complex domain.

% SECTION 4: shows the results of the optimization, including the search
% path layered over the cost fucntion and the corresponding learning curve.

clearvars; close all; clc;

%% Visualize physical system

%%%%%%%%%% Autonomous functions for physical caluclations %%%%%%%%%%

% Physical parameters of system
load parameterData.mat

% Complex numbers corresponding to coordinates in the foil-domian
z = @(z1,z2) z1 + 1i*z2;
x = @(z1,z2) 0.5*(z(z1,z2) + sign(z1).*(z(z1,z2).^2 - 4*k^2).^0.5);

% Location of image vortex (models vortex interaction with foil)
xim = @(z1v,z2v) x0 + r0^2./conj(x(z1v,z2v)-x0);

% Strength of circulation needed to meet Kutta condition
Gamma_0 = @(z1v,z2v) Gamma_v*((abs(x(z1v,z2v)-x0)^2 - r0^2)/(abs(x(z1v,z2v)-k)^2)) + ...
                  4*pi*u0*((delta-k)*sin(alpha) - beta*cos(alpha));

% Complex flow potential function
F = @(z1,z2,z1v,z2v) u0*exp(-1i*alpha)*x(z1,z2) + u0*exp(1i*alpha)*r0^2./(x(z1,z2)-x0) + ...
                     Gamma_0(z1v,z2v)*log(x(z1,z2)-x0)/(2*pi*1i) - ...
                     Gamma_v*log(x(z1,z2)-xim(z1v,z2v))/(2*pi*1i) + ...
                     Gamma_v*log(x(z1,z2)-x(z1v,z2v))/(2*pi*1i);

% Complex velocity field around cylinder
V = @(z1,z2,z1v,z2v) u0*exp(1i*alpha) - u0*exp(-1i*alpha)*r0^2./conj((x(z1,z2)-x0).^2) + ...
                     1i*Gamma_0(z1v,z2v)./(2*pi*conj(x(z1,z2)-x0)) + ...
                     1i*Gamma_v./(2*pi*conj(x(z1,z2)-x(z1v,z2v))) - ...
                     1i*Gamma_v./(2*pi*conj(x(z1,z2)-xim(z1v,z2v)));

% Complex velocity field around foil
W = @(z1,z2,z1v,z2v) V(z1,z2,z1v,z2v).*conj((x(z1,z2).^2)./(x(z1,z2).^2  - k^2));

% Real-valued pressure around the foil (kPa)
P = @(z1,z2,z1v,z2v) -rho*abs(W(z1,z2,z1v,z2v)).^2/2000;

%%%%%%%%%% Flow velocity calculations %%%%%%%%%%

% Choice of vortex location for visualization
zv = 2 + 2i;

% Shape definition of cylinder and foil
angle = linspace(0,2*pi,1000);
circle = r0*exp(1i*angle) + x0;
foil = circle + k^2./circle;

% Coordinate spacing
nGridX = 800; nGridY = 400;     % grid density
quiverSpaceX = 1:20:nGridX;     % quiver density in x-direction
quiverSpaceY = 1:10:nGridY;     % quiver density in y-direction

% Create meshgrid of coordinates
z1 = linspace(-8*k0,8*k0,nGridX); 
z2 = linspace(-4*k0,4*k0,nGridY); 
[Z1,Z2] = meshgrid(z1,z2); 
Z = z(Z1,Z2); X = x(Z1,Z2); 

% Eliminate points inside the cylinder/foil
Z(abs(X - x0) < r0) = NaN; X(abs(X - x0) < r0) = NaN;
X1 = real(X); X2 = imag(X);
Z1 = real(Z); Z2 = imag(Z);

% Compute complex potential function (to show velocity fields)
Phi = F(Z1,Z2,real(zv),imag(zv));

% Compute velocity fields directly
xDot = V(Z1,Z2,real(zv),imag(zv)); xDot = xDot(quiverSpaceY,quiverSpaceX);
zDot = W(Z1,Z2,real(zv),imag(zv)); zDot = zDot(quiverSpaceY,quiverSpaceX);
xDot = xDot./abs(xDot); zDot = zDot./abs(zDot); % make quiver arrows unit length

% Compute pressure field around foil
pressure = P(Z1,Z2,real(zv),imag(zv)); % pressure in kPa
vortexEps = 0.2; pressure(abs(Z - zv) < vortexEps) = NaN;

% Compute forces on the foil (optional)
surfacePressure = P(real(foil),imag(foil),real(zv),imag(zv));
drag = -cumtrapz(imag(foil),surfacePressure); drag = drag(end);
lift = cumtrapz(real(foil),surfacePressure); lift = lift(end);

% Location of vortex in cylinder domain
xv = x(real(zv),imag(zv));

%%%%%%%%%% Flow field plot %%%%%%%%%%

fig1 = figure(1); % define figure and size
set(fig1,'units','normalized','outerposition',[0.1 0.1 0.8 0.7])
contourf(Z1,Z2,imag(Phi),20); hold on;
fill(real(foil),imag(foil),'k')
scatter(real(zv),imag(zv),200,'or','Linewidth',2)
quiver(Z1(quiverSpaceY,quiverSpaceX),Z2(quiverSpaceY,quiverSpaceX),real(zDot),imag(zDot),'k')
axis equal; hold off;
legend('Flow Streamlines','','Vortex Location','Unit Direction of Flow',...
       'interpreter','latex','fontsize',16)
xlabel('In-Stream Distance (m)','interpreter','latex','fontsize',20)
ylabel('Cross-Stream Distance (m)','interpreter','latex','fontsize',20)
title('Flow Field Around Foil (m/s)','interpreter','latex','fontsize',20)

%%%%%%%%%% Pressure plot %%%%%%%%%%

fig2 = figure(2); % define figure and size
set(fig2,'units','normalized','outerposition',[0.1 0.1 0.8 0.7])
contourf(Z1,Z2,pressure,50); hold on;
fill(real(foil),imag(foil),'k')
scatter(real(zv),imag(zv),100,'or','Linewidth',2)
colorbar; axis equal; hold off;
xlabel('In-Stream Distance (m)','interpreter','latex','fontsize',20)
ylabel('Cross-Stream Distance (m)','interpreter','latex','fontsize',20)
title('Pressure Field Around Foil (kPa)','interpreter','latex','fontsize',20)

%% Cost functions (timely! avoid running multiple times)

clearvars; clc;

% Load physical parameters
load parameterData.mat

% Foil interior with boundary
angle = linspace(0,2*pi,1000);
circle = r0*exp(1i*angle) + x0;
foil = circle + k^2./circle;

% Coordinate spacing
nGridX = 400; nGridY = 200; % grid density

% Create meshgrid of coordinates
zv1 = linspace(-8*k0,8*k0,nGridX); 
zv2 = linspace(-4*k0,4*k0,nGridY); 
[ZV1,ZV2] = meshgrid(zv1,zv2); 
ZV = ZV1 + 1i*ZV2; XV = 0.5*(ZV + sign(ZV1).*(ZV.^2 - 4*k^2).^0.5);

% Eliminate points inside the cylinder/foil
% exemptArea = inpolygon(ZV1,ZV2,real(offsetFoil),imag(offsetFoil));
exemptArea = inpolygon(ZV1,ZV2,real(foil),imag(foil));
ZV(exemptArea) = NaN;
ZV1 = real(ZV); ZV2 = imag(ZV);

% Create storage for data
dragSpace = zeros(nGridY,nGridX);
penaltySpace = zeros(nGridY,nGridX);

% Iterate over space to calculate drag and penalty
for ii = 1:nGridY
    disp([num2str(100*ii/nGridY),'% complete'])
    for jj = 1:nGridX
        if ~isnan(ZV(ii,jj))
            dragSpace(ii,jj) = dragFunction([ZV1(ii,jj);ZV2(ii,jj)]);
            penaltySpace(ii,jj) = penaltyFunction([ZV1(ii,jj);ZV2(ii,jj)]);
        else
            dragSpace(ii,jj) = nan;
            penaltySpace(ii,jj) = nan;
        end
    end
end

% Sum two values for the objective function
objectiveSpace = dragSpace + penaltySpace;

% Save values to avoid re-running this section
save costFunctions.mat ZV1 ZV2 dragSpace penaltySpace objectiveSpace

%% Objective function components
% Note: you must run main.m before running this section!

clearvars; clc;

% Load data
load parameterData.mat
load outputData.mat
load costFunctions.mat

% Shape definition of cylinder and foil
angle = linspace(0,2*pi,1000);
circle = r0*exp(1i*angle) + x0;
foil = circle + k^2./circle;

% Truncate largest values in each drag to make visualization easier
dragTrunc = 0.8;
dragSpaceVis = dragSpace; 
dragSpaceVis(dragSpace > dragTrunc) = dragTrunc; 
dragSpaceVis(dragSpace < -dragTrunc) = -dragTrunc;

% Truncate largest values in penalty array to make visualization easier
penaltyTrunc = 1;
penaltySpaceVis = penaltySpace; 
penaltySpaceVis(penaltySpace > penaltyTrunc) = penaltyTrunc; 
penaltySpaceVis(penaltySpace < -penaltyTrunc) = -penaltyTrunc;

% Plot the drag function and penalty function seperately
fig3 = figure(3); % define figure and size
set(fig3,'units','normalized','outerposition',[0.1 0.1 1 0.6])

subplot(1,2,1)
contourf(ZV1,ZV2,dragSpaceVis,25); hold on; %,'edgecolor','none'); hold on;
fill(real(foil),imag(foil),'k')
hold off; axis equal; colorbar
xlabel('In-Stream Distance (m)','interpreter','latex','fontsize',18)
ylabel('Cross-Stream Distance (m)','interpreter','latex','fontsize',18)
title('Drag at a Given Vortex Position (kN)','interpreter','latex','fontsize',18)

subplot(1,2,2)
contourf(ZV1,ZV2,penaltySpaceVis,50,'edgecolor','none'); hold on; %); hold on;
fill(real(foil),imag(foil),'k')
hold off; axis equal; colorbar
xlabel('In-Stream Distance (m)','interpreter','latex','fontsize',18)
ylabel('Cross-Stream Distance (m)','interpreter','latex','fontsize',18)
title('Penalty at a Given Vortex Position (kN)','interpreter','latex','fontsize',18)

%% Results

clearvars; clc;

% Load data
load parameterData.mat
load outputData.mat
load costFunctions.mat

% Shape definition of cylinder and foil
angle = linspace(0,2*pi,1000);
circle = r0*exp(1i*angle) + x0;
foil = circle + k^2./circle;

% Truncate largest values in objective array to make visualization easier
objectiveTrunc = 1;
objectiveSpaceVis = objectiveSpace; 
objectiveSpaceVis(objectiveSpace > objectiveTrunc) = objectiveTrunc; 
objectiveSpaceVis(objectiveSpace < -objectiveTrunc) = -objectiveTrunc;

% Plot the objective function with the minimization path
fig4 = figure(4); % define figure and size
set(fig4,'units','normalized','outerposition',[0.1 0.1 1 1])
contourf(ZV1,ZV2,objectiveSpaceVis,50); hold on; %,'edgecolor','none'); hold on;
fill(real(foil),imag(foil),'k'); hold on;
plot(iterX(1,:),iterX(2,:),'-r','linewidth',2)
scatter(iterX(1,1),iterX(2,1),250,'or','linewidth',2)
scatter(iterX(1,end),iterX(2,end),250,'xr','linewidth',2)
hold off; axis equal; colorbar
xlabel('In-Stream Distance (m)','interpreter','latex','fontsize',24)
ylabel('Cross-Stream Distance (m)','interpreter','latex','fontsize',24)
legend('','Foil Outline','Search Path','Initial Search Position',...
       'Local Minimum Returned','interpreter','latex','fontsize',18)
title('Objective Function at a Given Vortex Position (kN)','interpreter','latex','fontsize',24)

% Plot the learning curve
fig5 = figure(5);
set(fig5,'units','normalized','outerposition',[0.1 0.1 0.5 0.6])
plot(0:length(iterF)-1,iterF,'-b','linewidth',1.5); hold on;
scatter(0:length(iterF)-1,iterF,75,'ob','filled'); hold on;
grid on; grid minor;
xlabel('Iteration Number','interpreter','latex','fontsize',18)
ylabel('Cost Function Evaluation','interpreter','latex','fontsize',18)
title('Quasi-Newton Algorithm Learning Curve','interpreter','latex','fontsize',18)
