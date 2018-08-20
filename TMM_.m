% Homework #3, Problem #1
% EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% This MATLAB program implements the transfer matrix method.
% INITIALIZE MATLAB
close all;
clc;
clear all;
% UNITS
degrees = pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
lam0 = 2.7; %free space wavelength
theta = 57 * degrees; %elevation angle
phi = 23 * degrees; %azimuthal angle
pte = 1/sqrt(2); %amplitude of TE polarization
ptm = 1i/sqrt(2); %amplitude of TM polarization
% EXTERNAL MATERIALS
ur1 = 1.2; %permeability in the reflection region
er1 = 1.4; %permittivity in the reflection region
ur2 = 1.6; %permeability in the transmission region
er2 = 1.8; %permittivity in the transmission region
% DEFINE LAYERS
UR = [ 1 3 ]; %array of permeabilities in each layer
ER = [ 2 1 ]; %array of permittivities in each layer
L = [ .25 .5 ]; %array of the thickness of each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE I AND 0 FOR MODES
I = eye(2);
O = zeros(2);
% LHI SIMPLE CASE
ninc1 = sqrt(er1 * ur1);
kx = 1.0006;% sin(theta) * cos(phi);
ky = 0.4247; %sin(theta) * sin(phi);
Qg = [  kx * ky   1 + ky^2;
        -(1+kx^2)  -kx * ky;];
Vg = -1i * Qg;
Leig = 1i;
% INITIALIZE S GLOBAL
S11g = O; S12g = I; S22g = O; S21g= I; 
% LOOP THROUGH ALL LAYERS
N_layer = size(ER,2);
for i = 1:N_layer
    % CALCULATE PARAMETERS FOR LAYER I
    kz(i) = sqrt(ER(i) * UR(i) - kx^2 - ky^2);
    Qi = 1 / UR(i) * [kx * ky   UR(i) * ER(i) - kx^2;
                        ky^2 - UR(i) * ER(i)  -kx * ky;];
    OMi = 1i * kz(i) * eye(size(Qi));
    Vi = Qi / OMi;
    % CALCULATE SCATTERING MATRIX FOR LAYER I
    Ai = eye(size(Vi)) + inv(Vi) * Vg;
    Bi = eye(size(Vi)) - inv(Vi) * Vg;
    Xi = expm(OMi * ko * L(i));
    D = Ai - Xi * Bi \ Ai * Xi * Bi;
    S11 = inv(D) * (Xi * Bi \ Ai * Xi * Ai - Bi);
    S22 = S11;
    S12 = inv(D) * Xi * (Ai - Bi \ Ai * Bi);
    S21 = S12;
    % UPDATE GLOABAL SCATTERING MATRIX
    Dg = S12g \ (I - S11 * S22g);
    Fg = S21 \ (I - S22g * S11);
    S11g = S11g + Dg * S11 * S21g;
    S12g = Dg *  S12;
    S21g = Fg * S21g;
    S22g = S22 + F * S22g * S12;
    
end