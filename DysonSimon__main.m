% 3H2 + N2 <=> 2NH3
% Dyson-Simon's kinetic equation
clear all; close all; clc

p      = 200;         % pressure in PFR      [atm]
n0N2   = 1;           % inlet mol flow of N2 [mol/hod]
n0H2   = 3;           % inlet mol flow of H2 [mol/hod]
n0NH3  = 0.01;        % inlet mol flow o NH3 [mol/hod]
R      = 8.314;       % molar gas constant   [J/mol/K]
Tin    = 350+273.15;  % temp                 [K]
nyN2   = -1.0;        % stoich coef of N2    [-]
nyH2   = -3.0;        % stoich coef of H2    [-]
nyNH3  = +2.0;        % stoich coef of NH3   [-]

pars = [p n0N2 n0H2 n0NH3 R Tin nyN2 nyH2 nyNH3];

yini = [n0N2 n0H2 n0NH3 Tin];

[VVou,yyou] = ode45(@(V,y) DysonSimon_actWITHOUT(V,pars,y),[0 0.00000000004],yini);
yyou(:,4) = yyou(:,4)-273.15;

[VVDS,yyDS] = ode45(@(V,y) DysonSimon_actDS(V,pars,y),[0 0.000000000004], yini);
yyDS(:,4) = yyDS(:,4)-273.15; % [°C]


[VVPR,yyPR] = ode45(@(V,y) DysonSimon_actPR(V,pars,y),[0 0.000000000004], yini);
yyPR(:,4) = yyPR(:,4)-273.15; % [°C]


[VVUNI,yyUNI] = ode45(@(V,y) DysonSimon_actUNI(V,pars,y),[0 0.000000000004], yini);
yyUNI(:,4) = yyUNI(:,4)-273.15;

figure;
plot(VVou,yyou(:,1),VVou,yyou(:,2),VVou,yyou(:,3)); box on; grid on; legend('N2','H2','NH3');
figure;
plot(VVou,yyou(:,4)); legend('Temp');

figure;
plot(VVDS, yyDS(:,1), '-b' ,VVDS, yyDS(:,2), '-r' ,VVDS, yyDS(:,3), '-k' ,...
     VVPR, yyPR(:,1), '--b',VVPR, yyPR(:,2), '--r',VVPR, yyPR(:,3), '--k',...
     VVUNI,yyUNI(:,1),'-.b',VVUNI,yyUNI(:,2),'-.r',VVUNI,yyUNI(:,3),'-.k'); box on; grid on;
legend('nN2DS','nH2DS','nNH3DS',...
       'nN2PR','nH2PR','nNH3PR',...
       'nN2UN','nH2UN','nNH3UN');
xlabel('Volume [m3]'); ylabel('Mole Flow [mol/hod]');

figure;
plot(VVDS, yyDS(:,4), '-',...
     VVPR, yyPR(:,4), '--',...
     VVUNI,yyUNI(:,4),'-.');
legend('TempDS','TempPR','TempUN');
box on; grid on; xlabel('Volume [m3]'); ylabel('Temperature [°C]');

figure;
plot(VVPR,yyPR(:,1),'--b',VVPR,yyPR(:,2),'--r',VVPR,yyPR(:,3),'--k',...
     VVou,yyou(:,1),'-.b',VVou,yyou(:,2),'-.r',VVou,yyou(:,3),'-.k');
legend('N2PR','H2PR','NH3PR','N2withoutact','H2withoutact','NH3withoutact'); close all;

A = [(n0N2-yyou(end,1))/n0N2  0                       0;
     (n0N2-yyDS(end,1))/n0N2 (n0N2-yyPR(end,1))/n0N2 (n0N2-yyUNI(end,1))/n0N2]; disp(A);
 
B = [max(yyou(:,4)) 0              0;
     max(yyDS(:,4)) max(yyPR(:,4)) max(yyUNI(:,4))]; disp(B);

