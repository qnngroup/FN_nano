clear;
clc;
close all;

%%Plot Setup%%
format compact;
format short g;
more on;

FontSize = 20;    %12
LineWidth = 2;  %1

set(0, 'DefaultTextFontSize', FontSize);
set(0, 'DefaultAxesFontSize', FontSize);
set(0, 'DefaultAxesFontName', 'Verdana');
set(0, 'DefaultAxesLineWidth', LineWidth);
set(0, 'DefaultAxesTickLength', [0.02 0.025]);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultLineLineWidth', LineWidth);
set(0, 'DefaultLineMarkerSize', 5);
set(0, 'DefaultPatchLineWidth', LineWidth);
set(0, 'DefaultFigureColor', [1 1 1]);
set(0, 'DefaultFigurePaperPosition',[2.65 4.3 3.2 2.4]);
set(0, 'DefaultFigurePosition', [800 500 360 270])
% 360x260 pixels on screen keeps same aspect ratio as 3.2x2.4 printed
set(0, 'DefaultFigurePosition', [300 0 960 960])
set(0, 'DefaultFigureDockControls', 'off')
%%Plot Setup%%

%Load functions from master directory:
addpath('../../');

%Load other experimental data sets
Ns = 200;
Fexp = linspace(1/1,1/0.01,Ns); %Field [V/nm]
Rexp = [7.3 10.3 6.5];    %nm
Aeffexp = [1 1.5e-6 10.2]; %nm^2
WFexp = [4.35 4.05 4.5];  %eV
betaexp = [0.07 0.0165 1];    %nm^-1,F = beta*V
Ntipexp = [25e3 1.36e6 1];    %number of tips
data1 = load('Spindt2010.csv');
V1 = data1(:,1);
I1 = data1(:,2) * 1e9;  %[A] -> [nA]
Fdata1 = V1 * betaexp(1);  %[V/nm]
data2 = load('Guerrera2012.csv');
V2 = data2(:,1);
I2 = data2(:,2) * 1e9;  %[A] -> [nA]
Fdata2 = V2 * betaexp(2);  %[V/nm]
data3 = load('Cabrera2013.csv');    %Already in F-N plot format!
x3 = data3(:,1);
y3 = data3(:,2);
Iexp = zeros(length(Rexp),length(Fexp));
for i = 1:length(Rexp)
    Iexp(i,:) = Jnano(WFexp(i),Fexp,Rexp(i)) * 1e-18 * 1e9 * Aeffexp(i) * Ntipexp(i); %[A/m^2] -> [nA/nm^2]
end

%Load Marco diodes' data
fn1 = 'nvt72';
data1 = load([fn1 '.mat']);
Invt72 = data1.I_E * 1e9;
Vnvt72 = data1.V_C;
Gnvt72 = 25;    %Gap1 [nm]
% FEnvt72 = 5;    %Field enhancement 2
FEnvt72 = 10;    %Field enhancement 1
betanvt72 = FEnvt72/Gnvt72; %F = beta * V = FE * V / G, beta = FE / G [nm^-1]
Fnvt72 = betanvt72 * Vnvt72;    %[V/nm]

%Fit Marco diodes' data
Rnvt72 = 20;    %[nm]
Aeffnvt72 = 1000; %[nm^2]
% WFnvt72 = 2.9;    %[eV]
WFnvt72 = 4.6;    %[eV]
Ntipnvt72 = 1;    
Ifitnvt72 = Jnano(WFnvt72,Fexp,Rnvt72) * 1e-18 * 1e9 * Aeffnvt72 * Ntipnvt72; %[A/m^2] -> [nA/nm^2]


figure(1)
corder = get(gca,'colororder');
scatter(1./Fnvt72,log(Invt72./Fnvt72.^2));
hold on;
plot(1./Fexp,log(Ifitnvt72./Fexp.^2),'color',corder(1,:));
scatter(1./Fdata1,log(I1./Fdata1.^2),'r');
xlim([0 1]);
ylim([-10 15]);
xlabel('1/F (nm/V)');ylabel('log(I/F^2) (I in nA, F in V/nm)');
plot(1./Fexp,log(Iexp(1,:)./Fexp.^2),'r');
scatter(1./Fdata2,log(I2./Fdata2.^2),'b');
plot(1./Fexp,log(Iexp(2,:)./Fexp.^2),'b');
scatter(x3,y3,'k');
plot(1./Fexp,log(Iexp(3,:)./Fexp.^2),'k');
legend({['NVT72,FE=' num2str(FEnvt72)], ['Wf=' num2str(WFnvt72) ' eV,R=' num2str(Rnvt72) ...
    ' nm,\beta=' num2str(betanvt72) ' nm^-^1,A_e_f_f=' num2str(Aeffnvt72) ' nm^2'], ...
    'Spindt et al. 2010, 25\times10^3 tips','Wf=4.35 eV,R=7.3 nm,\beta=0.07 nm^-^1,A_e_f_f=1 nm^2', ...
    'Guerrera et al. 2012, 1.36\times10^6 tips','Wf=4.05 eV,R=10.3 nm,\beta=0.0165 nm^-^1,A_e_f_f=1.5e-6 nm^2', ...
    'Cabrera et al. 2013, 1 tip','Wf=4.5 eV,R=6.5 nm,\beta=unknown,A_e_f_f=10.2 nm^2'},...
    'Location','Southoutside','NumColumns',1);

% Child-Langmuir law, I ~ V^3/2
N_CL = 200;
V_CL = logspace(-1,2,N_CL);
K_CL = 40;   % C-L law coefficient; I = K * V^3/2
I_CL = K_CL * V_CL.^1.5;

figure(2)
loglog(Vnvt72,Invt72,'o');
xlabel('V (V)');ylabel('I (nA)');
xlim([10^-1 10^1.5]);
hold on;
loglog(V_CL,I_CL,'k--');
legend({['NVT72'],...
    ['Child-Langmuir law, I=KV^1^.^5, K=' num2str(K_CL)]},...
    'Location','Southoutside','NumColumns',1);

% Schottky emission, log(I) ~ V^1/2
figure(3)
plot(sqrt(Vnvt72),log(Invt72),'o');
xlabel('V^1^/^2 (V in V)');ylabel('ln(I) (I in nA)');
hold on;
legend({['NVT72']},...
    'Location','Southoutside','NumColumns',1);