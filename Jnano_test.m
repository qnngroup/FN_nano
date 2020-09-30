clear;
clc;

%%Plot Setup%%
format compact;
format short g;
more on;

FontSize = 15;    %12
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
set(0, 'DefaultFigurePosition', [300 100 960 540])
set(0, 'DefaultFigureDockControls', 'off')
%%Plot Setup%%

Ns = 200;   %Sampling points
% F = linspace(0.1,2,Ns);    % field [V/nm]
% WF = 5.1;  % work function [eV]
F = linspace(1/0.35,1/0.1,Ns);  %Field [V/nm]
WF = 4.5;   %Workfunction [eV]
R = [5 10 20];  %Tip ROC [nm]

J = zeros(length(R),length(F));
for i = 1:length(R)
    J(i,:) = Jnano(WF,F,R(i));
end

J = J * 1e-18;  %[A/m^2] -> [A/nm^2]

figure(1)
plot(1./F,log(J(1,:)./F.^2),'k');
ylim([-40 -15]);
xlabel('1/F (nm/V)');ylabel('log(J/F^2) (J in A/nm^2, F in V/nm)');
hold on;
plot(1./F,log(J(2,:)./F.^2),'g');
plot(1./F,log(J(3,:)./F.^2),'r');
legend('R = 5 nm','R = 10 nm','R = 20 nm');

F2 = linspace(1/0.27,1/0.1,Ns); %Field [V/nm]
R2 = [7.3 10.3 6.5 4.2];    %nm
Aeff = [1 1.5e-6 10.2 6.3]; %nm^2
WF2 = [4.35 4.05 4.5 4.5];  %eV
% beta = [0.074 0.017 0.99 0.105];    %nm^-1,F = beta*V
beta = [0.07 0.0165 1 0.105];    %nm^-1,F = beta*V
Ntip = [25e3 1.36e6 1 1];    %number of tips

I = zeros(length(R2),length(F2));
for i = 1:length(R2)
    I(i,:) = Jnano(WF2(i),F2,R2(i)) * 1e-18 * 1e9 * Aeff(i) * Ntip(i); %[A/m^2] -> [nA/nm^2]
end

figure(2)
plot(1./F2,log(I(1,:)./F2.^2),'r');
xlim([0.1 0.27]);ylim([-10 15]);
xlabel('1/F (nm/V)');ylabel('log(I/F^2) (I in nA, F in V/nm)');
hold on;
plot(1./F2,log(I(2,:)./F2.^2),'b');
plot(1./F2,log(I(3,:)./F2.^2),'k');
plot(1./F2,log(I(4,:)./F2.^2),'g');
lgd = legend('R = 7.3 nm, A_e_f_f = 1 nm^2','R = 10.3 nm, A_e_f_f = 1.5\times10^-^6 nm^2', ...
    'R = 6.5 nm, A_e_f_f = 10.2 nm^2','R = 4.2 nm, A_e_f_f = 6.3 nm^2');
lgd.FontSize = 10;

data1 = load('Spindt2010.csv');
V1 = data1(:,1);
I1 = data1(:,2) * 1e9;  %[A] -> [nA]
Fdata1 = V1 * beta(1);  %[V/nm]
data2 = load('Guerrera2012.csv');
V2 = data2(:,1);
I2 = data2(:,2) * 1e9;  %[A] -> [nA]
Fdata2 = V2 * beta(2);  %[V/nm]
data3 = load('Cabrera2013.csv');    %Already in F-N plot format!
x3 = data3(:,1);
y3 = data3(:,2);

figure(3)
scatter(1./Fdata1,log(I1./Fdata1.^2),'r');
xlim([0.1 0.27]);ylim([-25 15]);
xlabel('1/F (nm/V)');ylabel('log(I/F^2) (I in nA, F in V/nm)');
hold on;
plot(1./F2,log(I(1,:)./F2.^2),'r');
scatter(1./Fdata2,log(I2./Fdata2.^2),'b');
plot(1./F2,log(I(2,:)./F2.^2),'b');
scatter(x3,y3,'k');
plot(1./F2,log(I(3,:)./F2.^2),'k');
legend('Spindt et al. 2010, 25e3 tips','Wf=4.35 eV,R=7.3 nm,\beta=0.07 nm^-^1,A_e_f_f=1 nm^2', ...
    'Guerrera et al. 2012, 1.36e6 tips','Wf=4.05 eV,R=10.3 nm,\beta=0.0165 nm^-^1,A_e_f_f=1.5e-6 nm^2', ...
    'Cabrera et al. 2013, 1 tip','Wf=4.5 eV,R=6.5 nm,\beta=unknown,A_e_f_f=10.2 nm^2',...
    'Location','Southwest');
