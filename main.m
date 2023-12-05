%%  Optimization of vortex position

clearvars; close all; clc;

% Physical parameters of system
load parameterData.mat

% Define interior penalty objective function
objectiveFun = @(z) penaltyFunction(z) + dragFunction(z);

% Sample initial guesses for the uncambered foil (beta = 0.0)
z0 = [-3;1];  % finds local minimum at the leading edge
% z0 = [4;1];   % finds local minimum at the trailing edge
% z0 = [-4;-1]; % finds otherwise hard to find local minimum

% Sample initial guesses for the cambered foil (beta = 0.2)
% z0 = [-2;-2];       % finds local minimum near the leading edge
% z0 = [1.65;-0.47];  % finds local minimum at the trailing edge
% z0 = [1;-0.5];      % finds otherwise hard to find local minimum

% Choose optimization algorithm parameters and solve 
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
          'FiniteDifferenceType','central','OutputFcn',@outputFunction);
[zopt,fval,exitflag,output] = fminunc(objectiveFun,z0,options);

% Set a global variable to store outputs at each iteration
global outputGlobal
iterX = horzcat(outputGlobal.x);
iterF = horzcat(outputGlobal.eval);

% Save optimization data as a .mat file
save outputData.mat iterX iterF
