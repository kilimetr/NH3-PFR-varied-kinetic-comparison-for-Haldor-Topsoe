function [dydz] = DysonSimon_actDS(V,pars,yvec)
% 3H2 + N2 <=> 2NH3
% Dyson-Simon's kinetic equation

p      = pars(1);  % pressure in PFR      [atm]
n0N2   = pars(2);  % inlet mol flow of N2 [mol/s]
n0H2   = pars(3);  % inlet mol flow of H2 [mol/s]
n0NH3  = pars(4);  % inlet mol flow o NH3 [mol/s]
R      = pars(5);  % molar gas constant   [J/mol/K]
Tin    = pars(6);  % temp                 [K]
nyN2   = pars(7);  % stoich coef of N2    [-]
nyH2   = pars(8);  % stoich coef of H2    [-]
nyNH3  = pars(9);  % stoich coef of NH3   [-]

nN2  = yvec(1);
nH2  = yvec(2);
nNH3 = yvec(3);
T    = yvec(4);

k  = 8.849 * 10^14 * exp(-40765 / (R*T));                % reaction velocity constant
KA = 10^(-2.691122 * log10(T) - 5.519265*10^(-5)*T + ...
     1.848863*10^(-7)*T^2 + 2001.6/T + 2.6899);          % equilibrium velocity constant

nTOT = nN2 + nH2 + nNH3;

yN2  = nN2  / nTOT;
yH2  = nH2  / nTOT;
yNH3 = nNH3 / nTOT;

phiN2  = 0.93431737 + 0.3101804*10^(-3)*T + 0.295896*10^(-3)*p - ...
         0.2707279*10^(-6)*T^2 + 0.4775207*10^(-6)*p^2;
phiH2  = exp(exp(-3.8402*T^0.125 + 0.541)*p - exp(-0.1263*T^0.5 - 15.98)*p^2 ...
         + 300*exp(-0.011901*T - 5.941)*exp(-p/300) - 1);
phiNH3 = 0.1438996 + 0.2028538*10^(-2)*T - 0.4487672*10^(-3)*p - ...
         0.1142945*10^(-5)*T^2 + 0.2761216*10^(-6)*p^2;

fN2  = yN2  * phiN2  * p * T/Tin;
fH2  = yH2  * phiH2  * p * T/Tin;
fNH3 = yNH3 * phiNH3 * p * T/Tin;

cpH2   = 28.660 +  1.170e-3 *T -  0.920e-6  *T^2;  % [J/mol/K]
cpN2   = 27.820 +  4.180e-3 *T;
cpNH3  = 27.31  + (2.383e-2)*T + (0.1707e-4)*T^2 + (-1.185e-8)*T^3;

drH = p*(0.5426 + 840.609/T + 4.5973*10^8/T^3) - 5.34685*T - 2.525*10^(-4)/T^2 + ...
               1.69167*10^(-6)*T^3 - 9157.09; % [Kcal/kmol]
drH = drH * 4.184; % [J/mol]

r = k*(KA^2 * fN2 * (fH2^3/fNH3^2)^0.5 - (fNH3^2/fH2^3)^0.5);

dnN2dV  = nyN2  * r;
dnH2dV  = nyH2  * r;
dnNH3dV = nyNH3 * r;

dTdV = -r*drH;
dTdV = dTdV / (nH2*cpH2 + nN2*cpN2 + nNH3*cpNH3);

dydz = [dnN2dV dnH2dV dnNH3dV dTdV]';

end
