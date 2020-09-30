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

%Load Marco diodes' data
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
G1 = 12;    %Gap1 [nm]
G2 = 27;    %Gap2 [nm]
FE1 = 40;    %Field enhancement 1
FE2 = 40;    %Field enhancement 2
beta1 = FE1/G1; %F = beta * V = FE * V / G, beta = FE / G [nm^-1]
beta2 = FE2/G2;
F1pos = beta1 * V1pos;  %[V/nm]
F1neg = beta1 * V1neg;  %[V/nm]
F2pos = beta2 * V2pos;  %[V/nm]
F2neg = beta2 * V2neg;  %[V/nm]

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

figure(1)
scatter(1./F1pos,log(I1pos./F1pos.^2));
hold on;
scatter(1./F1neg,log(I1neg./F1neg.^2));
scatter(1./F2pos,log(I2pos./F2pos.^2));
scatter(1./F2neg,log(I2neg./F2neg.^2));
scatter(1./Fdata1,log(I1./Fdata1.^2),'r');
xlim([0 2]);
ylim([-10 15]);
xlabel('1/F (nm/V)');ylabel('log(I/F^2) (I in nA, F in V/nm)');
plot(1./Fexp,log(Iexp(1,:)./Fexp.^2),'r');
scatter(1./Fdata2,log(I2./Fdata2.^2),'b');
plot(1./Fexp,log(Iexp(2,:)./Fexp.^2),'b');
scatter(x3,y3,'k');
plot(1./Fexp,log(Iexp(3,:)./Fexp.^2),'k');
legend({['NVT45-10-02-m04 positive,FE=' num2str(FE1)],['NVT45-10-02-m04 negative,FE=' num2str(FE1)], ...
    ['NVT45-02-01-m03 positive,FE=' num2str(FE2)],['NVT45-02-01-m03 negative,FE=' num2str(FE2)], ...
    'Spindt et al. 2010, 25e3 tips','Wf=4.35 eV,R=7.3 nm,\beta=0.07 nm^-^1,A_e_f_f=1 nm^2', ...
    'Guerrera et al. 2012, 1.36e6 tips','Wf=4.05 eV,R=10.3 nm,\beta=0.0165 nm^-^1,A_e_f_f=1.5e-6 nm^2', ...
    'Cabrera et al. 2013, 1 tip','Wf=4.5 eV,R=6.5 nm,\beta=unknown,A_e_f_f=10.2 nm^2'},...
    'Location','Southoutside','NumColumns',1);

% Child-Langmuir law, I ~ V^3/2
N_CL = 200;
V_CL = logspace(-1,1,N_CL);
K_CL = 1;   % C-L law coefficient; I = K * V^3/2
I_CL = K_CL * V_CL.^1.5;

figure(2)
loglog(V1pos,I1pos,'o');
xlabel('V (V)');ylabel('I (nA)');
xlim([10^-1 10^1]);
hold on;
loglog(V1neg,I1neg,'o');
loglog(V2pos,I2pos,'o');
loglog(V2neg,I2neg,'o');
loglog(V_CL,I_CL,'k--');
legend({['NVT45-10-02-m04 positive'],['NVT45-10-02-m04 negative'], ...
    ['NVT45-02-01-m03 positive'],['NVT45-02-01-m03 negative'], ...
    ['Child-Langmuir law, I=KV^1^.^5, K=' num2str(K_CL)]},...
    'Location','Southoutside','NumColumns',1);

% Schottky emission, log(I) ~ V^1/2
figure(3)
plot(sqrt(V1pos),log(I1pos),'o');
xlabel('V^1^/^2 (V in V)');ylabel('ln(I) (I in nA)');
hold on;
plot(sqrt(V1neg),log(I1neg),'o');
plot(sqrt(V2pos),log(I2pos),'o');
plot(sqrt(V2neg),log(I2neg ),'o');
legend({['NVT45-10-02-m04 positive'],['NVT45-10-02-m04 negative'], ...
    ['NVT45-02-01-m03 positive'],['NVT45-02-01-m03 negative'], ...
    },...
    'Location','Southoutside','NumColumns',1);