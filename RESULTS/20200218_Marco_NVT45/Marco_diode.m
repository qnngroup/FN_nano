clear;
clc;

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
set(0, 'DefaultFigurePosition', [300 100 960 540])
set(0, 'DefaultFigureDockControls', 'off')
%%Plot Setup%%

%Load functions from master directory:
addpath('../../');

fn1 = 'NVT45_10_02_m04';
fn2 = 'NVT45_02_01_m03';

data1 = load([fn1 '.mat']);
data2 = load([fn2 '.mat']);

I1 = data1.I_C * 1e9;
V1 = data1.V_C;
I2 = data2.I_C * 1e9;
V2 = data2.V_C;
V1pos = V1(V1>=0);
I1pos = I1(V1>=0);
V1neg = -V1(V1<0);
I1neg = -I1(V1<0);
V2pos = V2(V2>=0);
I2pos = I2(V2>=0);
V2neg = -V2(V2<0);
I2neg = -I2(V2<0);

figure(1)
plot(V1,I1);
xlabel('V (V)');ylabel('I (nA)');
hold on;
plot(V1pos,I1pos);
plot(V1neg,I1neg);
legend('raw','positive V','negative V','Location','Southeast');

figure(2)
plot(V2,I2);
xlabel('V (V)');ylabel('I (nA)');
hold on;
plot(V2pos,I2pos);
plot(V2neg,I2neg);
legend('raw','positive V','negative V','Location','Southeast');

G1 = 12;    %Gap1 [nm]
G2 = 27;    %Gap2 [nm]
Ns = 500;
F = linspace(1/1,1/0.01,Ns); %Field [V/nm]
R1 = 1000;    %nm
R2 = 1000;
Aeff1 = 5e-4; %nm^2
Aeff2 = 5;
% Aeff1 = 5e-4; %nm^2
% Aeff2 = 1000;
WF1 = 5.1;  %eV
WF2 = 5.1;
FE1 = 20;    %nm^-1,F = beta*V/G
FE2 = 25;
% FE1 = 20;    %nm^-1,F = beta*V/G
% FE2 = 10;

Jfit1 = Jnano(WF1,F,R1) * 1e-18;  %[A/m^2] -> [A/nm^2]
Ifit1 = Jfit1 * Aeff1 * 1e9;  %[A] -> [nA]
Vfit1 = F * G1 / FE1;

Jfit2 = Jnano(WF2,F,R2) * 1e-18;  %[A/m^2] -> [A/nm^2]
Ifit2 = Jfit2 * Aeff2 * 1e9;  %[A] -> [nA]
Vfit2 = F * G2 / FE2;

figure(3)
scatter(1./V1pos,log(I1pos./V1pos.^2));
xlim([0.12 0.15]);
ylim([-5 -2]);
xlabel('1/V_C (1/V)');ylabel('log(I_C/V_C^2) (I in nA, V_C in V)');
hold on;
scatter(1./V1neg,log(I1neg./V1neg.^2));
plot(1./Vfit1,log(Ifit1./Vfit1.^2));
legend('NVT45-10-02-m04 positive','NVT45-10-02-m04 negative', ...
    ['Wf = ' num2str(WF1) ' eV,' 'ROC = ' num2str(R1) ' nm,' ...
    'A_e_f_f = ' num2str(Aeff1) ' nm^2,' ' FE = ' num2str(FE1)], ...
    'Location','Southwest');

figure(4)
scatter(1./V2pos,log(I2pos./V2pos.^2));
% xlim([0.16 0.21]);
% ylim([-5 -2]);
xlabel('1/V_C (1/V)');ylabel('log(I_C/V_C^2) (I in nA, V_C in V)');
hold on;
scatter(1./V2neg,log(I2neg./V2neg.^2));
plot(1./Vfit2,log(Ifit2./Vfit2.^2));
legend('NVT45-02-01-m03 positive','NVT45-02-01-m03 negative', ...
    ['Wf = ' num2str(WF2) ' eV,' 'ROC = ' num2str(R2) ' nm,' ...
    'A_e_f_f = ' num2str(Aeff2) ' nm^2,' ' FE = ' num2str(FE2)], ...
    'Location','Southwest');
