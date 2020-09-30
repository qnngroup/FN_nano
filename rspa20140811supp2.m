function [J]=GenFNEqFun(F,W,R)
%This function implements the general model for field emission. The inputs
%are arrays with the same dimensions. 
%The file SpecialFunctions.mat must be in the current folder in order to
%load the special functions as computed algebraically.
%The input arguments are F: local electric field in V/nm, W: work function
%in eV, R: radius of curvature at the apex of the emitter in nm.
%If one of the input variables is an array, the other input variables must
%be either arrays of the same dimensions or scalars.

yy=1.2*sqrt(F)./W;%the standard variable y


M=csvread('SpecialFunctions.csv');%load the special functions data file
y=M(:,1);
v=M(:,2);
t=M(:,3);
w=M(:,4);
ps=M(:,5);

vv=interp1(y,v,yy);%interpolating special function
ww=interp1(y,w,yy);
tt=interp1(y,t,yy);
ppsi=interp1(y,ps,yy);

G=6.831*(vv.*(W.^1.5)./F+ww.*(W.^2.5)./(R.*F.^2));%Gamow exponent

d=(sqrt(W)./F).*(tt+(W./(F.*R)).*ppsi);%prefactor

J=1.54e-6*(d.^(-2)).*exp(-G);%calculation of current density (in A/nm^2)

end