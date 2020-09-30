function J = Jnano(WF,F,R)
% Compute the current density with Fowler-Nordheim theory on a nanoscale
% tip.
% Input: WF (work function [eV]), F (electric field [V/nm]), R (tip ROC [nm])
% Output: J (current density [A/m^2])
% All units are in SI.

%FN constants
e_const = 1.602176462e-19;  % electron charge [C]
eps0_const = 8.854187817e-12;   % vacuum permittivity[F/m]
me = 9.10938188e-31;    % electron mass [kg]
h = 6.62607004e-34;     % Planck constant [m^2 kg s^-1]
aFN = e_const.^3/(8*pi*h);
bFN = 8*pi*sqrt(2*me)/(3*e_const*h);

%Unit conversion
WF = WF * e_const;  %[eV] -> [J]
R = R * 1e-9;   %[nm] -> [m]
F = F * 1e9;    %[V/nm] -> [V/m]

%FN corrections
B = e_const^3/(16*pi*eps0_const);   % [C^2 F^-1 m]
y = 2*sqrt(B.*F)./WF;
% vy = 1 - y.^2 - (y.^2.*log(y))/3;
% wy = 4/5 - 7*y.^2/40 - (y.^2.*log(y))/100;
% ty = 1 + y.^2/9 - (y.^2.*log(y))/11;
% psiy = 4/3 - y.^2/500 - (y.^2.*log(y))/15;
%correction functions in tabular form; cite:
%Kyritsakis, A. and Xanthakis, J.P., 2015. Derivation of a generalized
%Fowler–Nordheim equation for nanoscopic field-emitters. Proceedings of the
%Royal Society A: Mathematical, Physical and Engineering Sciences,
%471(2174), p.20140811.
taby = csvread('rspa20140811supp3.csv');     
vy = interp1(taby(:,1),taby(:,2),y);
ty = interp1(taby(:,1),taby(:,3),y);
wy = interp1(taby(:,1),taby(:,4),y);
psiy = interp1(taby(:,1),taby(:,5),y);

J = aFN.*WF.^(-1).*F.^2.*(ty+WF./(e_const.*F.*R).*psiy).^(-2) ...
    .*exp(-bFN.*WF.^(3/2)./F.*(vy+WF./(e_const.*F.*R).*wy));
end