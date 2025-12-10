% Spec Heat Capacity: 'Cp', 'C' (/1000)
% Density: 'Ro', 'D'
% Internal Energy: 'u', 'U' (/1000)
% Enthalpy: 'h', 'H' (/1000)
% Entropy: 's', 'S' (/1000)
% Viscosity: 'ViscDyn', 'V'
% Therm Conductivity: 'ThermCond', 'L' (/1000)

% H2O:    bad fit for ViscDyn liquid water <50°C
% R134a:  bad fit for Cp, ViscDyn

T = (10:10:400) + 273.15;
P = [0.1; 1; 10; 100] * 101325;
S = 5000:1000:10000;
fluid = 'H2O';

%{
resultA = zeros(size(T));    % size(T) = [1 21]
resultB = resultA;
for j=1:length(P)
    for i=1:length(T)           % length(T) = 21
    resultA(i) = ADVENS_Props('Cp', 'P', P(j),'T',T(i), fluid);
    resultB(i) = 1000*Props('C', 'P', P(j)/1000,'T',T(i), fluid);
    end
    figure
    plot(T, [resultA; resultB])
    legend('ADVENS','CoolProp')
end
%}

resultA = zeros(size(T));    % size(T) = [1 21]
resultB = resultA;
for j=1:length(S)
    for i=1:length(T)           % length(T) = 21
    resultA(i) = ADVENS_Props('Ro', 's', S(j),'T',T(i), fluid);
    resultB(i) = Props('D', 'S', S(j)/1000,'T',T(i), fluid);
    end
    figure
    plot(T, [resultA; resultB])
    legend('ADVENS','CoolProp')
end