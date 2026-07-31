%% ========================================================================
% Resonant Tunneling Through a Double Barrier Structure
% NEGF + Tight Binding Method
%
% Structure:
% Lead | Barrier(10 nm,25 eV) | Well(20 nm) | Barrier(10 nm,25 eV) | Lead
%
% Outputs:
%   1. Potential Profile
%   2. Transmission Spectrum T(E)
%   3. First Resonant State Density
%   4. Second Resonant State Density
% ========================================================================

clear;
clc;
close all;

%% ------------------------------------------------------------------------
% Physical Constants
%% ------------------------------------------------------------------------

hbar = 1.054571817e-34;
m0   = 9.1093837015e-31;
q    = 1.602176634e-19;

% hbar^2/(2m0) in eV·nm^2
C = hbar^2/(2*m0)/(q*1e-18);

%% ------------------------------------------------------------------------
% Device Parameters
%% ------------------------------------------------------------------------

a  = 0.25;      % lattice spacing (nm)

Lc = 10;        % lead length (nm)
Lb = 10;        % barrier width (nm)
Lw = 20;        % well width (nm)

Vb = 25;        % barrier height (eV)

m  = m0;

%% ------------------------------------------------------------------------
% Tight-Binding Hopping Parameter
%% ------------------------------------------------------------------------

t0 = C*(m0/m)/a^2;

fprintf('t0 = %.4f eV\n',t0);

%% ------------------------------------------------------------------------
% Potential Profile
%% ------------------------------------------------------------------------

segL = [Lc Lb Lw Lb Lc];
segV = [0  Vb  0  Vb  0];

U = [];
for k = 1:length(segL)
    U = [U segV(k)*ones(1,round(segL(k)/a))]; %#ok<AGROW>
end

N = length(U);
x = (0:N-1)*a;

%% ------------------------------------------------------------------------
% Hamiltonian
%% ------------------------------------------------------------------------

H = diag(U + 2*t0) ...
  + diag(-t0*ones(N-1,1),1) ...
  + diag(-t0*ones(N-1,1),-1);

I = eye(N);

%% ------------------------------------------------------------------------
% Transmission Spectrum
%% ------------------------------------------------------------------------

EmeV = linspace(0.2,6,5000);
EeV = EmeV/1000;

T = zeros(size(EeV));

for n = 1:length(EeV)
    E = EeV(n);

    ka = acos(1 - E/(2*t0));
    sigma = -t0*exp(1i*ka);

    Sigma1 = zeros(N);
    Sigma2 = zeros(N);

    Sigma1(1,1) = sigma;
    Sigma2(N,N) = sigma;

    Gamma1 = 1i*(Sigma1 - Sigma1');
    Gamma2 = 1i*(Sigma2 - Sigma2');

    G = (E*I - H - Sigma1 - Sigma2)\I;

    T(n) = real(trace(Gamma1*G*Gamma2*G'));
end

%% ------------------------------------------------------------------------
% Plot Transmission
%% ------------------------------------------------------------------------

figure('Color','w');
plot(EmeV,T,'b','LineWidth',2);

grid on;
box on;

xlabel('Energy (meV)');
ylabel('Transmittance T(E)');
title('Transmission Spectrum');

%% ------------------------------------------------------------------------
% Resonant States (isolated well estimate)
%% ------------------------------------------------------------------------

wellStart = round((Lc+Lb)/a)+1;
wellEnd   = round((Lc+Lb+Lw)/a);

Hw = H(wellStart:wellEnd,wellStart:wellEnd);

[V,D] = eig(Hw); %#ok<ASGLU>
[Eres,idx] = sort(real(diag(D))); %#ok<ASGLU>

E1 = Eres(1);
E2 = Eres(2);

fprintf('\n');
fprintf('First Resonance  E1 = %.4f meV\n',E1*1000);
fprintf('Second Resonance E2 = %.4f meV\n',E2*1000);

%% ------------------------------------------------------------------------
% Ideal Resonant-State Density (for clean visualization)
%% ------------------------------------------------------------------------

xwell_local = linspace(0,Lw,1000);

rho1_local = sin(pi*xwell_local/Lw).^2;
rho2_local = sin(2*pi*xwell_local/Lw).^2;

rho1_local = rho1_local/max(rho1_local);
rho2_local = rho2_local/max(rho2_local);

figure('Color','w','Position',[100 100 900 650]);

subplot(2,1,1)
plot(xwell_local,rho1_local,'r','LineWidth',3);

grid on;
box on;

xlabel('Position x (nm)');
ylabel('|\psi_1(x)|^2');

title(sprintf('First Resonant State  (E_1 = %.3f meV)',E1*1000));

xlim([0 Lw]);
ylim([0 1.1]);

subplot(2,1,2)
plot(xwell_local,rho2_local,'b','LineWidth',3);

grid on;
box on;

xlabel('Position x (nm)');
ylabel('|\psi_2(x)|^2');

title(sprintf('Second Resonant State (E_2 = %.3f meV)',E2*1000));

xlim([0 Lw]);
ylim([0 1.1]);

%% ------------------------------------------------------------------------
% Potential Profile
%% ------------------------------------------------------------------------

figure('Color','w');
plot(x,U,'k','LineWidth',2);

grid on;
box on;

xlabel('Position (nm)');
ylabel('Potential (eV)');

title('Double Barrier Potential Profile');

%% ------------------------------------------------------------------------
% Barrier Profile + Resonant States
%% ------------------------------------------------------------------------

figure('Color','w','Position',[100 100 1000 500])

plot(x,U,'k','LineWidth',3)
hold on

yline(E1,'r--','LineWidth',2)
yline(E2,'b--','LineWidth',2)

xwell_abs = linspace(Lc+Lb, Lc+Lb+Lw, 1000);

psi1 = sin(pi*(xwell_abs-(Lc+Lb))/Lw);
psi2 = sin(2*pi*(xwell_abs-(Lc+Lb))/Lw);

psi1 = psi1/max(abs(psi1));
psi2 = psi2/max(abs(psi2));

scale = 0.6; % visual scaling

plot(xwell_abs, E1 + scale*psi1,'r','LineWidth',2)
plot(xwell_abs, E2 + scale*psi2,'b','LineWidth',2)

xline(Lc,'k:')
xline(Lc+Lb,'k:')
xline(Lc+Lb+Lw,'k:')
xline(Lc+2*Lb+Lw,'k:')

grid on
box on

xlabel('Position x (nm)')
ylabel('Energy (eV)')

title('Double-Barrier Potential Profile and Resonant States')

legend('Potential Profile', ...
    sprintf('E_1 = %.4f meV',E1*1000), ...
    sprintf('E_2 = %.4f meV',E2*1000), ...
    'Location','northwest')

%% ------------------------------------------------------------------------
% Potential Profile + Resonant States + Electron Density
%% ------------------------------------------------------------------------

figure('Color','w','Position',[100 100 1000 650]);

yyaxis left
plot(x,U,'k','LineWidth',3);
hold on;

xlabel('Position x (nm)');
ylabel('Potential Energy (eV)');
ylim([0 Vb*1.1]);

yline(E1,'r--','LineWidth',2);
yline(E2,'b--','LineWidth',2);

legend('Potential Profile', ...
    sprintf('E_1 = %.3f meV',E1*1000), ...
    sprintf('E_2 = %.3f meV',E2*1000), ...
    'Location','northwest');

yyaxis right

rho1_abs = sin(pi*(xwell_abs-(Lc+Lb))/Lw).^2;
rho2_abs = sin(2*pi*(xwell_abs-(Lc+Lb))/Lw).^2;

rho1_abs = rho1_abs/max(rho1_abs);
rho2_abs = rho2_abs/max(rho2_abs);

plot(xwell_abs,rho1_abs,'r','LineWidth',3);
plot(xwell_abs,rho2_abs,'b','LineWidth',3);

ylabel('Normalized Electron Density');
ylim([0 1.1]);

xline(Lc,'k:');
xline(Lc+Lb,'k:');
xline(Lc+Lb+Lw,'k:');
xline(Lc+2*Lb+Lw,'k:');

grid on;
box on;

title('Double Barrier Potential Profile and Resonant-State Electron Densities');
