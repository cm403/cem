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
% .theta elevation angle of incidence (radians)
% .phi azimuthal angle of incidence (radians)
% .ate amplitude of TE polarization
% .atm amplitude of TM polarization


e0 = 8.854187817e-12 ;
u0 = pi * 4e-7;
c0 = 1/sqrt(u0*e0);
Z0 = sqrt(u0/e0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
src.lam0 = 1; %free space wavelength 
% k0 = 2 * pi /lam0
src.theta = 0; %elevation angle
src.phi = 0; %azimuthal angle
src.pte = 1;%1/sqrt(2); %amplitude of TE polarization
src.ptm = 0;%1/sqrt(2); %amplitude of TM polarization

% EXTERNAL MATERIALS
dev.ur1 = 1 * ( Z0 / Z0)^2; %permeability in the reflection region
dev.er1 = 1 / ( Z0 / Z0)^2; %permittivity in the reflection region
dev.ur2 = 1 * ( Z0 / Z0)^2; %permeability in the transmission region
dev.er2 = 1 / ( Z0 / Z0)^2; %permittivity in the transmission region
% DEFINE LAYERS
dev.UR = [ 1 1 1 1]; %array of permeabilities in each layer
dev.ER = [ 1 80 2 1]; %array of permittivities in each layer
% dev.ER = [ 1 / ( sqrt(50*100) / Z0)^2]; %array of permittivities in each layer

% % EXTERNAL MATERIALS
% dev.ur1 = 1 * ( 50 / 377)^2; %permeability in the reflection region
% dev.er1 = 1; %permittivity in the reflection region
% dev.ur2 = 1 * ( 100 / 377)^2; %permeability in the transmission region
% dev.er2 = 1; %permittivity in the transmission region
% % DEFINE LAYERS
% dev.UR = [ 1 * ( 70.7 / 377)^2;]; %array of permeabilities in each layer
% dev.ER = [ 1]; %array of permittivities in each layer

Zg = sqrt(dev.ur1/dev.er1 * u0/e0);
ZL = sqrt(dev.ur2/dev.er2 * u0/e0);
Z1 = 70.7;

disp(sprintf('\nZg = %.f', Zg));
disp(sprintf('Z1 = %.f', Z1));
disp(sprintf('ZL = %.f', ZL));

f = (0:.01:20) * 1e9;

Pref = nan(size(f));
Ptrn = nan(size(f));
S11 = nan(size(f));
S12 = nan(size(f));
S21 = nan(size(f));
S22 = nan(size(f));

% dev.L = c0 / 10e9 * .25 * ( 1 / sqrt(dev.UR * dev.ER) );
dev.L = [ 2 2 2 2] ;

for i = 50%1:size(f,2)
    
    src.lam0 = c0 / f(i);
%     dat = tmm1d(dev,src);
    dat = tmm1d_fields(dev,src);
    Pref(i) = dat.REF;
    Ptrn(i) = dat.TRN;
    S11(i) = dat.S11(1);
    S21(i) = dat.S21(1);
    S12(i) = dat.S12(1);
    S22(i) = dat.S22(1);
end

figure(1)
clf
subplot(2,1,1)
hold on
% plot(f/1e9, 20*log10(abs(S11)))
% plot(f/1e9, 20*log10(abs(S21)))
plot(f/1e9, 10*log10(Pref), 'DisplayName', 'Ref')
plot(f/1e9, 10*log10(Ptrn), 'DisplayName', 'Trn')
ylim([-60 0])
xlabel('Frequency [GHz]')
ylabel('S11 [dB20]')
grid on
legend('show','location','southeast')

subplot(2,1,2)
[~, h1] = smithchart(S11);
h1.SubLineType %= ':'
hold on
% plot(S11,'*')
% plot(S22,'*')


Psi_minus = dat.Psi_minus;
Psi_plus = dat.Psi_plus;
Psi = dat.Psi;
zdist = [];
flds = [];
figure(3)
clf
for j = 0:-.1:-10*pi
    zdist = [];
    flds = [];
    flds_minus = [];
    flds_plus = [];
    for i = 1:size(Psi,3)
        lami =  c0 * 1/sqrt(dev.UR(i)*dev.ER(i))/f(50);
        txt = sprintf('\\lambda = %.2f',lami);

    %     subplot(size(Psi,3), 1, i)
        if i == 1
            z_off = 0;
        else
            z_off = z_off + dev.L(i-1);
        end
        flds_minus = [flds_minus real(exp(1i*j)*Psi_minus(2,:,i))];
        flds_plus = [flds_plus real(exp(1i*j)*Psi_plus(2,:,i))];
        flds = [flds real(exp(1i*j)*Psi(2,:,i))];
        
        zdist = [zdist linspace(z_off,dev.L(i) + z_off,100)];
        plot([z_off z_off], [-6 6],'k--','Linewidth',3)
        hold on
        x1 = pi;
        y1 = sin(pi);
        txtmed = sprintf('Medium %.f\n \\epsilon_r: %.1f\n\\mu_r: %.1f',i, dev.ER(i), dev.UR(i));
        text(z_off + dev.L(i)/2,5,txtmed,'HorizontalAlignment','center')
    %     plot(zdist,real(exp(1i*j)*Psi(2:3,:,i))')
    %     hold on
    %     title(txt)
    %     ylim([-6 6])
    end
    subplot(3,1,1)
    plot(zdist,flds,'k')
    ylim([-6 6])
    subplot(3,1,1)
    plot(zdist,flds_plus)
    ylim([-6 6])
    subplot(3,1,1)
    plot(zdist,flds_minus)
    hold off
    ylim([-6 6])
    pause(0.1)
end
