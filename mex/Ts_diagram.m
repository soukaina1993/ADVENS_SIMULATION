function plt = Ts_diagram(fluid)
hold off
fluid = convertStringsToChars(fluid);
data = readtable('DataFluids.txt','PreserveVariableNames',true,'ReadRowNames',true);
if strcmp(fluid, 'H2O')     % CoolProp used
    Tmin = 273.16;
    Tmax = 647.096;
else
    Tmin = 250;
    Tmax = data{fluid, 'Tc'};
end

x = [0 1];
T = linspace(Tmin, Tmax);
for j=1:length(T)
    for i=1:length(x)
        s_x(i,j) = ADVENS_Props('s', 'x', x(i), 'T', T(j), fluid);
    end
end
for i=1:length(x)
    plot(s_x(i,:), T, 'red')
    title('T,s-diagram')
    xlabel('s')
    ylabel('T')
    hold on
end