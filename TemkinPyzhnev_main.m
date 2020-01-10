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

[VV,yy] = ode45(@(V,y) TemkinPyzhnev(V,pars,y),[0 0.000000000003], yini);

yy(:,4) = yy(:,4)-273.15; % [°C]

figure(1);
subplot(2,1,1); plot(VV,yy(:,1)); hold on; box on; grid on; plot(VV,yy(:,2)); plot(VV,yy(:,3)); 
legend('nN2','nH2','nNH3'); xlabel('Volume [m3]'); ylabel('Mole Flow [mol/hod]');
subplot(2,1,2); plot(VV,yy(:,4)); box on; grid on; xlabel('Volume [m3]'); ylabel('Temp [°C]');

figure(2);
plot(VV,yy(:,4)); xlabel('Volume [m3]'); 
ylabel('Temperature [°C]'); box on; grid on;

disp((n0N2-yy(end,1))/n0N2);
