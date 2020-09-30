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

%FN constants
e_const = 1.602176462e-19;  % electron charge [C]
eps0_const = 8.854187817e-12;   % vacuum permittivity[F/m]
me = 9.10938188e-31;    % electron mass [kg]
h = 6.62607004e-34;     % Planck constant [m^2 kg s^-1]
aFN = e_const.^3/(8*pi*h);
bFN = 8*pi*sqrt(2*me)/(3*e_const*h);

%User defined parameters 
% N = 1000;   % sampling point
% F = linspace(0.1e9,2e9,N);    % field [V/m]
WF = 5.1;  % work function [eV]
WF = WF*e_const;  % work function [J]
% R1 = 50e-9;    % nanotip ROC [m]
% R2 = 20e-9;    % nanotip ROC [m]
a = 50; %half foci distance [nm]
V0 = 10; %voltage across the gap [V]

%cartesian coordinate
Xmin = -200;
Xmax = 200;
Ymin = -200;
Ymax = 200;
Zmin = -200;
Zmax = 200;
Nx = 100;
Ny = 100;
Nz = 100;
x = linspace(Xmin,Xmax,Nx);
y = linspace(Ymin,Ymax,Ny);
z = linspace(Zmin,Zmax,Nz);

%prolate spheroidal coordinate
% xi = 1/(2*a)*(sqrt(x.^2+y.^2+(z+a).^2)+sqrt(x.^2+y.^2+(z-a).^2));
% eta = 1/(2*a)*(sqrt(x.^2+y.^2+(z+a).^2)-sqrt(x.^2+y.^2+(z-a).^2));
% phi = atan(y./x);

%constant eta
ximax = 200;
xieps = 1e-9;   %numerical zero
ximin = 1+xieps;
Nxi = 100;
% xi = linspace(ximin,ximax,Nxi);
xi = 1 + logspace(log10(ximin-1),log10(ximax-1),Nxi);   %logspace better resolve the apex shape
Nphi = 200;
phi = linspace(0,2*pi,Nphi);
eta0 = 0.5;
[Xi,Phi] = meshgrid(xi,phi);
x0 = a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*cos(Phi);
y0 = a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*sin(Phi);
z0 = a.*Xi.*eta0;

%calculate ROC
% R = sqrt(abs(a.^2.*eta0.^2.*Xi.^3.*sqrt(Xi.^2-1).*(1./eta0.^2-1./Xi.^2).^(3/2)));
% min(R(:))

figure()
s = surf(x0,y0,z0);
xlim([-500 500]);ylim([-500 500]);zlim([0 200]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
% s.EdgeColor = 'none';
% colorbar
caxis([0 200]);

% partial derivatives of parameterized surface
rx = cat(3,a.*sqrt(1-eta0.^2).*cos(Phi).*Xi./sqrt(Xi.^2-1), ...
    a.*sqrt(1-eta0.^2).*sin(Phi).*Xi./sqrt(Xi.^2-1), ...
    a.*eta0.*ones(size(Xi)));
rp = cat(3,-a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*sin(Phi), ...
    a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*cos(Phi), ...
    zeros(size(Xi)));
rx_x_rp = cross(rx,rp);
rx_x_rp_norm = vecnorm(rx_x_rp,2,3);    %2-norm of the 3-rd dimension
n = rx_x_rp./rx_x_rp_norm;
rxx = cat(3,-a.*sqrt(1-eta0.^2).*cos(Phi).*(Xi.^2-1).^(-3/2), ...
    -a.*sqrt(1-eta0.^2).*sin(Phi).*(Xi.^2-1).^(-3/2), ...
    zeros(size(Xi)));
rxp = cat(3,-a.*sqrt(1-eta0.^2).*sin(Phi).*Xi./sqrt(Xi.^2-1), ...
    a.*sqrt(1-eta0.^2).*cos(Phi).*Xi./sqrt(Xi.^2-1), ...
    zeros(size(Xi)));
rpp = cat(3,-a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*cos(Phi), ...
    -a.*sqrt(Xi.^2-1).*sqrt(1-eta0.^2).*sin(Phi), ...
    zeros(size(Xi)));

% fundamental forms of the surface
E = dot(rx,rx,3);
F = dot(rx,rp,3);
G = dot(rp,rp,3);
L = dot(rxx,n,3);
M = dot(rxp,n,3);
N = dot(rpp,n,3);

% surface curvatures
K = (L.*N - M.^2)./(E.*G - F.^2);
H = (E.*N - 2.*F.*M + G.*L)./2./(E.*G - F.^2);

% [K,H,Pmax,Pmin] = surfature(x0,y0,z0);
RK = sqrt(1./K);
RH = 1./H;
deltaR = RK - RH;
RK_SI = RK / 1e9;   %nm to m

% tip electric field
F1 = -2*V0./a./sqrt(Xi.^2 - eta0.^2)./sqrt(1 - eta0.^2)./log((1 - eta0)./(1 + eta0));   %V/nm
F1_SI = F1 * 1e9;   %V/nm to V/m

%FN corrections
B = e_const^3/(16*pi*eps0_const);   % [C^2 F^-1 m]
y = 2*sqrt(B.*F1_SI)./WF;
vy = 1 - y.^2 - (y.^2.*log(y))/3;
wy = 4/5 - 7*y.^2/40 - (y.^2.*log(y))/100;
ty = 1 + y.^2/9 - (y.^2.*log(y))/11;
phiy = 4/3 - y.^2/500 - (y.^2.*log(y))/15;

% nano-tip FN
JnanoK = aFN.*WF.^(-1).*F1_SI.^2.*(ty+WF./(e_const.*F1_SI.*RK_SI).*phiy).^(-2).*exp(-bFN.*WF.^(3/2)./F1_SI.*(vy+WF./(e_const.*F1_SI.*RK_SI).*wy));

figure()
surf(x0,y0,RK);
xlim([-50 50]);ylim([-50 50]);zlim([0 300]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('ROC (Gaussian) (nm)');
caxis([0 300]);

figure()
surf(x0,y0,RH);
xlim([-50 50]);ylim([-50 50]);zlim([0 300]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('ROC (mean) (nm)');
caxis([0 300]);

figure()
surf(x0,y0,deltaR);
xlim([-50 50]);ylim([-50 50]);zlim([0 100]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('ROC difference (nm)');
caxis([0 100]);

figure()
yyaxis left;
plot(x0(1,:),RK(1,:));
xlabel('x (nm)');ylabel('ROC (nm)');
xlim([0 50]);ylim([-10 300]);
hold on;
plot(x0(1,:),RH(1,:));
yyaxis right;
plot(x0(1,:),deltaR(1,:));
ylabel('ROC difference (nm)');
ylim([-1 30]);
legend('Gaussian ROC','mean ROC','ROC difference','Location','NorthWest');

figure()
surf(x0,y0,F1,'EdgeColor','none');
xlim([-500 500]);ylim([-500 500]);zlim([0 0.5]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('tip electric field (V/nm)');
caxis([0 0.5]);

figure()
surf(x0,y0,JnanoK,'EdgeColor','none');
xlim([-50 50]);ylim([-50 50]);%zlim([0 0.06]);
xlabel('x (nm)');ylabel('y (nm)');zlabel('emission current (A/m^2)');
%caxis([0 0.06]);