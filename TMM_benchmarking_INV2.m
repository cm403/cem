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
% k0 = 2 * pi /lam0
theta = 0 * degrees; %elevation angle
phi = 23 * degrees; %azimuthal angle
pte = 1/sqrt(2); %amplitude of TE polarization
ptm = 1i/sqrt(2); %amplitude of TM polarization
% EXTERNAL MATERIALS
ur1 = 1.2; %permeability in the reflection region
er1 = 1.4; %permittivity in the reflection region
ur2 = 1.2; %permeability in the transmission region
er2 = 1.4; %permittivity in the transmission region
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
kxt = sqrt(ur1*er1) * sin(theta) * cos(phi);
kyt = sqrt(ur1*er1) * sin(theta) * sin(phi);
QG = [  kxt * kyt   1 + kyt^2;
        -(1+kxt^2)  -kxt * kyt;];
VG = -1i * QG;

% INITIALIZE S GLOBAL
S11G = O; S12G = I; S22G = O; S21G = I; 
% LOOP THROUGH ALL LAYERS
N_layer = size(ER,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BENCHMARKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k0 = 2*pi/lam0
kx = k0 * kxt
ky = k0 * kyt
QG
VG


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAYERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kzti = nan(size(ER));
i = 1;
for i = 1:N_layer
    % CALCULATE PARAMETERS FOR LAYER I
    kzti(i) = sqrt(ER(i) * UR(i) - kxt^2 - kyt^2);
    Qi = 1 / UR(i) * [kxt * kyt   UR(i) * ER(i) - kxt^2;
                        kyt^2 - UR(i) * ER(i)  -kxt * kyt;]
    OMEGAi = 1i * kzti(i) * eye(size(Qi))
    Vi = Qi / OMEGAi
    % CALCULATE SCATTERING MATRIX FOR LAYER I
    Ai = eye(size(Vi)) + Vi \ VG
    Bi = eye(size(Vi)) - Vi \ VG
%     [Eig_v, ~] = eig(Qi)
    Xi = expm(OMEGAi * k0 * L(i) * lam0)   
    D = Ai - Xi * (Bi / Ai) * Xi * Bi
    S11 = D \ (Xi * (Bi / Ai) * Xi * Ai - Bi)
    
    S12 = D \ (Xi * (Ai - (Bi / Ai) * Bi))
    S21 = S12
    S22 = S11
    
    % UPDATE GLOABAL SCATTERING MATRIX
%     Dg = S12g \ (I - S11 * S22g)
    DG = S12G / (I - S11 * S22G)
%     Fg = S21 \ (I - S22g * S11)
    FG = S21 / (I - S22G * S11)
    
    S11G = S11G + DG * S11 * S21G
    S12G = DG *  S12
    S21G = FG * S21G
    S22G = S22 + FG * S22G * S12
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTERNAL REGIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CALCULATE SCATTERING MATRIX FOR REFLECTION REGION
kztR = sqrt(er1 * ur1 - kxt^2 - kyt^2);
QR = 1 / ur1 * [kxt * kyt   ur1 * er1 - kxt^2;
                    kyt^2 - ur1 * er1  -kxt * kyt;]
OMEGAR = 1i * kztR * eye(size(QR))
VR = QR / OMEGAR

AR = eye(size(VR)) + VG \ VR
BR = eye(size(VR)) - VG \ VR


S11R = -AR \ BR 
S12R = 2 * inv(AR) 
S21R = 0.5 * ( AR - (BR / AR) * BR)
S22R = BR / AR

%% CALCULATE SCATTERING MATRIX FOR TRANSMISSION REGION
kztT = sqrt(er2 * ur2 - kxt^2 - kyt^2);
QT = 1 / ur2 * [kxt * kyt   ur2 * er2 - kxt^2;
                    kyt^2 - ur2 * er2  -kxt * kyt;]
OMEGAT = 1i * kztT * eye(size(QT))
VT = QT / OMEGAT
AT = eye(size(VT)) + VG \ VT
BT = eye(size(VT)) - VG \ VT


S11T = BT / AT
S12T = 0.5 * ( AT - (BT / AT) * BT)
S21T = 2 * inv(AT) 
S22T = -AT \ BT

%% CONNECT TO EXTERNAL REGIONS

% SG = SR ? SG
DR = S12R / (I - S11G * S22R)
FR = S21G / (I - S22R * S11G)

S22G = S22G + FR * S22R * S12G
S21G = FR * S21R
S12G = DR *  S12G
S11G = S11R + DR * S11G * S21R


% SG = SG ? ST
DG = S12G / (I - S11T * S22G)
FG = S21T / (I - S22G * S11T)

S11G = S11G + DG * S11T * S21G
S12G = DG *  S12T
S21G = FG * S21G
S22G = S22T + FG * S22G * S12T

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE SCATTERING PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if theta = 0
n = [0; 0; 1;];
kinc = k0 * ninc1 * ...
    [sin(theta) * cos(phi); sin(theta)*sin(phi); cos(theta)];
if theta == 0
    ate = [0 1 0];
else
    ate = cross(n, kinc) ./ norm(cross(n, kinc));
end
atm = cross(ate, kinc) ./ norm(cross(ate, kinc));
P  = pte * ate + ptm * atm;
Esrc = [P(1); P(2)]
Eref = [S11G * Esrc; nan]
Etrn = [S21G * Esrc; nan]

Eref(3) = -(kxt * Eref(1) + kyt * Eref(2)) / kztR
Etrn(3) = -(kxt * Etrn(1) + kyt * Etrn(2)) / kztT

REF = norm(Eref).^2
TRN = norm(Etrn).^2 * real(kztT * k0  / ur2) / real(kinc(3) / ur1)
