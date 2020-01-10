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

[VV,yy] = ode45(@(V,y) DysonSimon_actWITHOUT(V,pars,y),[0 0.000000000005], yini);

yy(:,4) = yy(:,4)-273.15; % [°C]

figure(1);
subplot(2,1,1); plot(VV,yy(:,1)); hold on; box on; grid on; plot(VV,yy(:,2)); plot(VV,yy(:,3)); 
legend('nN2','nH2','nNH3'); xlabel('Volume [m3]'); ylabel('Mole Flow [mol/hod]');
subplot(2,1,2); plot(VV,yy(:,4)); box on; grid on; xlabel('Volume [m3]'); ylabel('Temp [°C]');

figure;
plot(VV,yy(:,1),VV,yy(:,2),VV,yy(:,3)); box on; grid on;
legend('nN2','nH2','nNH3'); xlabel('Volume [m3]'); ylabel('Mole Flow [mol/hod]');

figure;
plot(VV,yy(:,4)); box on; grid on; xlabel('Volume [m3]'); ylabel('Temperature [°C]');

return;


% thermostuff

t = 20:5:600; tvar = t'; % [°C]
Tvar = tvar + 273.15;

for i=1:length(Tvar)
    k(i)  = 8.849 * 10^14 * exp(-40765 / (R*Tvar(i)));                      % reaction velocity constant
    KA(i) = 10^(-2.691122 * log10(Tvar(i)) - 5.519265*10^(-5)*Tvar(i) + ... % [mol i/m3/h]
            1.848863*10^(-7)*Tvar(i)^2 + 2001.6/Tvar(i) + 2.6899);          % equilibrium constant [1/atm]
end

figure(2); subplot(2,1,1); plot(tvar,k); xlabel('Temp [°C]');
ylabel('Forward Velocity Constant [mol "i"/m3/h*atm^(0.5)]'); box on; grid on;
subplot(2,1,2); plot(tvar,KA); xlabel('Temp [°C]'); ylabel('Equilibrium Velocity Constant [1/atm]'); 
box on; grid on;

for i=1:length(Tvar)
    T = Tvar(i);
    k1(i) = 1.78954*10^(4)*exp(-20800/(R*T));
    k2(i) = 2.5714*10^(16)*exp(-47400/(R*T));
    k22(i) = k(i)/KA(i);
end
    
figure;subplot(2,1,1); plot(tvar,k1,tvar,k); subplot(2,1,2); plot(tvar,k2,tvar,k22);


dH0H2  = 0;     % [J/mol]
dH0N2  = 0;
dH0NH3 = -45720;

dG0H2  = 0;     % [J/mol]
dG0N2  = 0;
dG0NH3 = -16400;


for i=1:length(Tvar)
    T = Tvar(i);
    
    cpH2(i)   = 28.660+ 1.170e-3*T- 0.920e-6*T^2;  % [J/mol/K]
    cpN2(i)   = 27.820+ 4.180e-3*T;
    cpNH3(i)  = 27.31+(2.383e-2)*T+(0.1707e-4)*T^2+(-1.185e-8)*T^3;
    
    drcp(i) = nyH2*cpH2(i) + nyN2*cpN2(i) + nyNH3*cpNH3(i);
    
    % activity coef by Dyson Simon
    phiDSN2(i)  = 0.93431737 + 0.3101804*10^(-3)*T + 0.295896*10^(-3)*p - ...
                0.2707279*10^(-6)*T^2 + 0.4775207*10^(-6)*p^2;
    
    phiDSH2(i)  = exp(exp(-3.8402*T^0.125 + 0.541)*p - exp(-0.1263*T^0.5 - 15.98)*p^2 ...
                + 300*exp(-0.011901*T - 5.941)*exp(-p/300) - 1);
    
    phiDSNH3(i) = 0.1438996 + 0.2028538*10^(-2)*T - 0.4487672*10^(-3)*p - ...
                0.1142945*10^(-5)*T^2 + 0.2761216*10^(-6)*p^2;
    
    t = Tvar(i) - 273.15;
    % activity coef taken from Aspen Plus described by UNIFAC at 100bar
    phiUNN2(i)  = -5*10^(-8)*t^2 + 4*10^(-5)*t + 1.0240;
    phiUNH2(i)  = +5*10^(-8)*t^2 - 9*10^(-5)*t + 1.0575;
    phiUNNH3(i) = -8*10^(-7)*t^2 + 0.0009   *t + 0.6958;
    % activity coef taken from Aspen Plus described by Peng-Robinson at
    % 100bar
    phiPRN2(i)  = +6*10^(-8)*t^2 + 3*10^(-5)*t + 1.0324;
    phiPRH2(i)  = +5*10^(-8)*t^2 - 8*10^(-5)*t + 1.0481;
    phiPRNH3(i) = -9*10^(-7)*t^2 + 0.0011   *t + 0.6716;
end


for i=2:length(Tvar)
    T = Tvar(i);
    
    dHH2(i)  = dH0H2  + cpH2(i) *(Tvar(i) - 293.15);
    dHN2(i)  = dH0N2  + cpN2(i) *(Tvar(i) - 293.15);
    dHNH3(i) = dH0NH3 + cpNH3(i)*(Tvar(i) - 293.15);
end


for i=1:length(Tvar)
    T = Tvar(i);
    drHDS(i) = p*(0.5426 + 840.609/T + 4.5973*10^8/T^3) - 5.34685*T - 2.525*10^(-4)/T^2 + ...
               1.69167*10^(-6)*T^3 - 9157.09; % [Kcal/kmol N3]
    drHDS(i) = drHDS(i) * 4.184;
end

drH = nyH2*dHH2 + nyN2*dHN2 + nyNH3*dHNH3;

figure(3);
subplot(2,2,1); plot(tvar,cpH2); hold on; box on; grid on; plot(tvar,cpN2); plot(tvar,cpNH3); 
legend('cpH2','cpN2','cpNH3'); xlabel('Temp [°C]'); ylabel('Heat Capacity [J/mol/K]');

subplot(2,2,2); plot(tvar,drcp); box on; grid on; xlabel('Temp [°C]'); 
ylabel('Reaction Heat Capacity [J/mol/K]');

subplot(2,2,3); plot(tvar,dHH2);  hold on; box on; grid on; plot(tvar,dHN2); plot(tvar,dHNH3); 
legend('dHH2','dHN2','dHNH3'); xlabel('Temp [°C]'); ylabel('Enthalpy [J/mol]');

subplot(2,2,4); plot(tvar,drH); box on; grid on; xlabel('Temp [°C]'); ylabel('Reaction Enthalpy [J/mol]');

figure(4);
plot(tvar,phiDSH2); hold on; box on; grid on; plot(tvar,phiDSN2); plot(tvar,phiDSNH3); 
legend('phiH2','phiN2','phiNH3'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients by Dyson-Simon');

figure(5);
plot(tvar,phiUNH2), hold on; box on; grid on; plot(tvar,phiUNN2), plot(tvar,phiUNNH3);
legend('phiH2','phiN2','phiNH3'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients Taken from Aspen Described by UNIFAC at 100 bar');

figure(6);
plot(tvar,phiPRH2), hold on; box on; grid on; plot(tvar,phiPRN2), plot(tvar,phiPRNH3);
legend('phiH2','phiN2','phiNH3'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients Taken from Aspen Described by Peng-Robinson at 100 bar');

figure(7);
plot(tvar,phiDSH2), hold on; box on; grid on; plot(tvar,phiUNH2); plot(tvar,phiPRH2);
legend('phiDS','phiUN','phiPR'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients for H2 Calculated Variously');

figure(8);
plot(tvar,phiDSN2), hold on; box on; grid on; plot(tvar,phiUNN2); plot(tvar,phiPRN2);
legend('phiDS','phiUN','phiPR'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients for N2 Calculated Variously');

figure(9);
plot(tvar,phiDSNH3), hold on; box on; grid on; plot(tvar,phiUNNH3); plot(tvar,phiPRNH3);
legend('phiDS','phiUN','phiPR'); xlabel('Temp [°C]'); ylabel('Activity coef [-]');
title('Activity Coefficients for NH3 Calculated Variously');

figure(10);
plot(tvar,drH); hold on; box on; grid on; plot(tvar,drHDS); legend('drH','drHDS');
xlabel('Temp [°C]'); ylabel('Reaction Enthalpy [kJ/kmol]'); ylim([min(drH) max(drHDS)]); 


figure; plot(tvar,phiDSH2,'-r', tvar,phiDSN2,'-b', tvar,phiDSNH3,'-k',...
             tvar,phiUNH2,'-.r',tvar,phiUNN2,'-.b',tvar,phiUNNH3,'-.k',...
             tvar,phiPRH2,'--r',tvar,phiPRN2,'--b',tvar,phiPRNH3,'--k'); box on; grid on;
legend('phiDSH2','phiDSN2','phiDSNH3',...
       'phiUNH2','phiUNN2','phiUNNH3',...
       'phiPRH2','phiPRN2','phiPRNH3');
title('Activity Coefficients - Variously Calculated'); xlabel('Temp [°C]'); ylabel('Activity [-]');

