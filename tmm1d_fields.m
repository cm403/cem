function DAT = tmm1d_fields(DEV,SRC)
    % TMM1D One-Dimensional Transfer Matrix Method
    %
    % DAT = tmm1d(DEV,SRC);
    %
    % INPUT ARGUMENTS
    % ================
    % DEV Device Parameters
    % .er1 relative permittivity in reflection region
    % .ur1 relative permeability in reflection region
    % .er2 relative permittivity in transmission region
    % .ur2 relative permeability in transmission region
    % .ER array containing permittivity of each layer
    % .UR array containing permeability of each layer
    % .L array containing thickness of each layer
    %
    % SRC Source Parameters
    % .lam0 free space wavelength
    % .theta elevation angle of incidence (degress)
    % .phi azimuthal angle of incidence (degrees)
    % .ate amplitude of TE polarization
    % .atm amplitude of TM polarization
    %
    % OUTPUT ARGUMENTS
    % ================
    % DAT Output Data
    % .REF Reflectance
    % .TRN Transmittance
    % .S11 GLOBAL (DEVICE)
    % .S12 GLOBAL (DEVICE)
    % .S21 GLOBAL (DEVICE)
    % .S22 GLOBAL (DEVICE)
    % This MATLAB program implements the transfer matrix method.
    
    % UNITS
    degrees = pi/180;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DEFINE SIMULATION PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOURCE PARAMETE0S
    lam0 = SRC.lam0; %free space wavelength 
    theta = SRC.theta * degrees; %elevation angle; convert from deg to rad
    phi = SRC.phi * degrees; %azimuthal angle; convert from deg to rad
    pte = SRC.pte; %amplitude of TE polarization
    ptm = SRC.ptm; %amplitude of TM polarization
    % EXTERNAL MATERIALS
    ur1 = DEV.ur1; %permeability in the reflection region
    er1 = DEV.er1; %permittivity in the reflection region
    ur2 = DEV.ur2; %permeability in the transmission region
    er2 = DEV.er2; %permittivity in the transmission region
    % DEFINE LAYERS
    UR = DEV.UR; %array of permeabilities in each layer
    ER = DEV.ER; %array of permittivities in each layer
    L = DEV.L; %array of the thickness of each layer
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
    WG = eye(2);
    VG = -1i * QG;
    k0 = 2 * pi / lam0;
    % INITIALIZE S GLOBAL
    S11G = O; S12G = I; S22G = O; S21G = I; 
    % CALCULATE # OF LAYERS FROM ER INPUT
    N_layer = size(ER,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LOOP THROUGH LAYERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    kzti = nan(size(ER));
    i = 1;
    for i = 1:N_layer
        % CALCULATE PARAMETERS FOR LAYER I
        kzti(i) = sqrt(ER(i) * UR(i) - kxt^2 - kyt^2);
        Qi = 1 / UR(i) * [kxt * kyt   UR(i) * ER(i) - kxt^2;
                            kyt^2 - UR(i) * ER(i)  -kxt * kyt;];
        OMEGAi = 1i * kzti(i) * eye(size(Qi));
        Vi = Qi / OMEGAi;
        % CALCULATE SCATTERING MATRIX FOR LAYER I
        Wi = eye(2);
        Ai = inv(Wi) * WG + inv(Vi) * VG;
        Bi = inv(Wi) * WG - inv(Vi) * VG;
        Xi = expm(OMEGAi * k0 * L(i) )   ;
        D = Ai - Xi * Bi / Ai * Xi * Bi;
        
        S11 = D \ (Xi * Bi / Ai * Xi * Ai - Bi);
        S12 = D \ Xi * (Ai - Bi / Ai * Bi);
        S21 = S12;
        S22 = S11;
        
        
        % UPDATE GLOABAL SCATTERING MATRIX
        DG = S12G * inv(I - S11 * S22G);
        FG = S21 * inv(I - S22G * S11);

        S11G = S11G + DG * S11 * S21G;
        S12G = DG *  S12;
        S21G = FG * S21G;
        S22G = S22 + FG * S22G * S12;

        S11G_f(:,:,i) = S11G + DG * S11 * S21G;
        S12G_f(:,:,i) = DG *  S12;
        S21G_f(:,:,i) = FG * S21G;
        S22G_f(:,:,i) = S22 + FG * S22G * S12;
        
        % SAVE INTERNAL FIELDS
        W(:,:,i) = eye(2);
        V(:,:,i) = Qi / OMEGAi;
        EIG(:,:,i) = OMEGAi;

    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXTERNAL REGIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % CALCULATE SCATTERING MATRIX FOR REFLECTION REGION
    kztR = sqrt(er1 * ur1 - kxt^2 - kyt^2);
    QR = 1 / ur1 * [kxt * kyt   ur1 * er1 - kxt^2;
                        kyt^2 - ur1 * er1  -kxt * kyt;];
    OMEGAR = 1i * kztR * eye(size(QR));
    VR = QR / OMEGAR;
    AR = eye(size(VR)) + inv(VG) * VR;
    BR = eye(size(VR)) - inv(VG) * VR;
    
    S11R = -AR \ BR ;
    S12R = 2 * inv(AR) ;
    S21R = 0.5 * ( AR - BR * inv(AR) * BR);
    S22R = BR / AR;
    
    % CALCULATE SCATTERING MATRIX FOR TRANSMISSION REGION
    kztT = sqrt(er2 * ur2 - kxt^2 - kyt^2);
    QT = 1 / ur2 * [kxt * kyt   ur2 * er2 - kxt^2;
                        kyt^2 - ur2 * er2  -kxt * kyt;];
    OMEGAT = 1i * kztT * eye(size(QT));
    VT = QT / OMEGAT;
    AT = eye(size(VT)) + inv(VG) * VT;
    BT = eye(size(VT)) - inv(VG) * VT;

    S11T = BT / AT;
    S12T = 0.5 * ( AT - BT * inv(AT) * BT);
    S21T = 2 * inv(AT); 
    S22T = -AT \ BT;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONNECT TO EXTERNAL REGIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SG = SR (X) SG
    DR = S12R * inv(I - S11G * S22R);
    FR = S21G * inv(I - S22R * S11G);

    S22G = S22G + FR * S22R * S12G;
    S21G = FR * S21R;
    S12G = DR *  S12G;
    S11G = S11R + DR * S11G * S21R;

    % SG = SG (X) ST
    DG = S12G * inv(I - S11T * S22G);
    FG = S21T * inv(I - S22G * S11T);

    S11G = S11G + DG * S11T * S21G;
    S12G = DG *  S12T;
    S21G = FG * S21G;
    S22G = S22T + FG * S22G * S12T;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE SCATTERING PROBLEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = [0; 0; 1;];
    kinc = k0 * ninc1 * ...
        [sin(theta) * cos(phi); sin(theta)*sin(phi); cos(theta)];
    
    % Check if theta = 0
    if theta == 0
        ate = [0 1 0];
    else
        ate = cross(n, kinc) ./ norm(cross(n, kinc));
    end
    atm = cross(ate, kinc) ./ norm(cross(ate, kinc));
    
    P  = pte * ate + ptm * atm;
    
    Esrc = [P(1); P(2)];
    Eref = [S11G * Esrc; nan];
    Etrn = [S21G * Esrc; nan];

    Eref(3) = -(kxt * Eref(1) + kyt * Eref(2)) / kztR;
    Etrn(3) = -(kxt * Etrn(1) + kyt * Etrn(2)) / kztT;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE MODE COEFF (BACKWARD PASS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cinc = [ P(1); P(2)];
    SG = [S11G S12G; S21G S22G];
    cref_trn = SG * [cinc; 0; 0;];
    cref = cref_trn(1:2);
    ctrn = cref_trn(3:4);
    for i = N_layer:-1:2
        z = linspace(0, L(i), 100);
        zp = k0 * z;
        c_minus(:,i) = S12G_f(:,:,i-1) \ (cref - S11G_f(:,:,i-1) * cinc);
        c_plus(:,i)  = S21G_f(:,:,i-1) * cinc + S22G_f(:,:,i-1) * c_minus(:,1);
        
        WV = [W(:,:,i), W(:,:,i);
              V(:,:,i), -V(:,:,i)];
        
        c = [c_plus(:,i); c_minus(:,i)];
        for zind = 1:size(zp,2)
            eEig = expm( [EIG(:,:,i) * zp(zind), zeros(2);...
                          zeros(2), -EIG(:,:,i) *  zp(zind)] );
            Psi(:,zind,i) = WV * eEig * c;
            Psi_plus(:,zind,i) = WV * eEig * [c_plus(:,i); zeros(size(c_minus(:,i)))];
            Psi_minus(:,zind,i) = WV * eEig * [zeros(size(c_plus(:,i))); c_minus(:,i)];
        end
        
    end
    
    for i = 1
        z = linspace(0, L(i), 100);
        zp = k0 * z;
        c_minus(:,i) = cref;
        c_plus(:,i)  = cinc;
        
        WV = [W(:,:,i), W(:,:,i);
              V(:,:,i), -V(:,:,i)];
        
        c = [c_plus(:,i); c_minus(:,i)];
        for zind = 1:size(zp,2)
            eEig = expm( [EIG(:,:,i) * zp(zind), zeros(2);...
                          zeros(2), -EIG(:,:,i) *  zp(zind)] );
            Psi(:,zind,i) = WV * eEig * c;
            Psi_plus(:,zind,i) = WV * eEig * [c_plus(:,i); zeros(size(c_minus(:,i)))];
            Psi_minus(:,zind,i) = WV * eEig * [zeros(size(c_plus(:,i))); c_minus(:,i)];
        end
    end
    
   
    DAT.REF = norm(Eref).^2;
    DAT.TRN = norm(Etrn).^2 * real(kztT * k0  / ur2) / real(kinc(3) / ur1);
    DAT.Eref = Eref;
    DAT.Etrn = Etrn * real(kztT * k0  / ur2) / real(kinc(3) / ur1);
    DAT.S11 = S11G;
    DAT.S12 = S12G;
    DAT.S21 = S21G;
    DAT.S22 = S22G;
    DAT.Psi = Psi;
    DAT.Psi_plus = Psi_plus;
    DAT.Psi_minus = Psi_minus;
    
    
end