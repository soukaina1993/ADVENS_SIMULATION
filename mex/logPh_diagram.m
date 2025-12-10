function plt = logPh_diagram(fluid)
hold off
fluid = convertStringsToChars(fluid);
data = readtable('DataFluids.txt','PreserveVariableNames',true,'ReadRowNames',true);
if strcmp(fluid, 'H2O')
    Pmin = 611.7;
    Pmax = 22064000;
elseif strcmp(fluid, 'R134a')
    Pmin = 100000;
    Pmax = data{fluid, 'Pc'};
else
    error("Not implemented")
end

x = [0 1];
P = logspace(log10(Pmin), log10(Pmax));
P(end) = Pmax;
for j=1:length(P)
    for i=1:length(x)
        h_x(i,j) = ADVENS_Props('h', 'x', x(i), 'P', P(j), fluid);
    end
end
for i=1:length(x)
    semilogy(h_x(i,:), P, 'red')
    title('P,h-diagram')
    xlabel('h')
    ylabel('P')
    hold on
end