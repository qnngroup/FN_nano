%This script fits the generalized FN equation to experimental data in order
%to extract the following parameters: enhancementa factor beta, radius of
%curvature R and prefactor sigma*Aefficient.

clear all;
close all;

W=4.5;%define the work function

%% give the experimental data straight from mat file
%There must be two vectors in the mat file: xexp and yexp. They contain the
%I-V data in the corresponding FN coordinates. xexp=1/V, yexp=log(I/V^2).
%Experimental data must be very smooth and accurate. Otherwise the method faces
%instabilities and it would be much better to fit a polynomial to the
%experimental data and give the data in polynomial form. In this case
%activate lines 26-30 in the next section.

load ExpData.mat;
xe=xexp;
ye=yexp;
%% give experimental data in polynomial form
% fit a polynomial to the data of the form yexp=p1*xexp^2+p2*xexp+p3
% activate this section if the data is noisy and you choose to fit a
% polynomial first

% [f1,g1]=fit(xexp,yexp,'poly2');%fit polynomial
% xstart=min(xexp); %give starting xexp
% xend=max(xexp); %give last xexp  
% xe=xstart:(xend-xstart)/50:xend;
% ye=f1(xe);%calculate values

%% main section

Ra=0.5;%extreme values for R
Rb=100;
for i=1:100
    R=0.5*(Ra+Rb);%bisection for R
    Fa=0.2;%extreme maximum fields
    Fb=20;
    for j=1:100
        Fmax=0.5*(Fa+Fb);%bisection for maximum field
        beta=Fmax*min(xe);%beta=Fmax/Vmax
        F=beta./xe;%local field vector corresponding to xexp vector
        I=GenFNEqFun(F,W,R);%current density
        yanal=log(I.*(xe.^2));%FNplot vertical variable
        sigma=max(yanal)-max(ye);%FN plot vertical offset to match
        yanal=yanal-sigma;%match plots at the left by appropriate prefactor
        if min(yanal)-min(ye)>0.001%
            Fb=Fmax;%bisection to the correct side, match plots at the right
        elseif min(ye)-min(yanal)>0.001%this way beta or Fmax is defined
            Fa=Fmax;
        else
            break;%if they match with 1% accuracy break the loop
        end
    end
    if mean(yanal)-mean(ye)>0.001
        Ra=R;%bisection at the correct side
    elseif mean(ye)-mean(yanal)>0.001
        Rb=R;%once left-right endpoints are matched, R is tuned to match the middle
    else
        break;%if the error at the mean is smaller than 0.001% break
    end
end
sigmaAeff=exp(-sigma);%correcting prefactor
%%   plotting the results
%get variables y for 1/F in x instead of 1/V
Fexp=beta./xexp;
Iexp=exp(yexp)./xexp.^2;
yexpF=log(Iexp./Fexp.^2);
Ianal=exp(yanal)./xe.^2;
yanF=log(Ianal./F.^2);

plot(1./Fexp,yexpF,'*');%plot experimental data with F in the free variable
hold on
plot(1./F,yanF,'r');%plot theoretical
xlabel('1/F (nm/V)');%label
ylabel('log(I/F^2) (I in A F in V/nm)');
figure;
plot(xexp,yexp,'*');%plot data with V in the free variable
hold on
plot(xe,yanal,'r');
xlabel('1/V (1/V)');
ylabel('log(I/V^2) (A/V^2)');
xlabel('1/V (1/V)');%label
ylabel('log(I/V^2) (A/V^2)');
disp(['The radius of curvature is R=' num2str(R) ' nm']);%show results of R, beta and sigmaAeff
disp(['The enhancement factor is beta=' num2str(beta) ' 1/nm']);
disp(['The pre-factor is sigmaAeff=' num2str(sigmaAeff) ' nm^2']);


