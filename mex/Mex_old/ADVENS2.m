%% ADVENS Thermodynamics : calcul de fluide
    % P [Pa]
    % T [K]
    
close all
clear all

cellArray{1}= "C2H60"; %water; R134a;
name='MEG-0%';
Ttp=230;
Tc=510;
T_test=350;
%P_test=300000;
choix_param= 1; % 1= (v,T); 2= (v,P); 3= (P,T); 4=(T,s); 5= (P,s); 6= (P,h), 7= Ajouter un fluide 
results= main(cellArray,[choix_param, 0.0019, Tc]);
Pc=results(18)-1;

choix_param= 1; % 1= (v,T); 2= (v,P); 3= (P,T); 4=(T,s); 5= (P,s); 6= (P,h), 7= Ajouter un fluide 
results= main(cellArray,[choix_param, 0.0019, T_test]);
P_test=results(18)-1;
%%
P0=50000;
x=0;


choix_param= 7; % 1= (v,T); 2= (v,P); 3= (P,T,x); 4=(T,s); 5= (P,s); 6= (P,h), P,x 7= Ajouter un fluide 
%P,T,v,u,h,s,cv,cp,x,sv,Ro,ViscDyn,ViscCin,ThermCond,Prandtl,vl,vg,Psat,ho; 
valeur_param1= 101325;
valeur_param2= 0;
results= main(cellArray,[choix_param, valeur_param1, valeur_param2]);


P=[P0:5000:Pc];
for i=1:length(P)
    results= main(cellArray,[choix_param, P(i),x]);
    h_0(i)=results(5);
    s_0(i)=results(6);
    T_0(i)=results(2);
    %h_0CP(i)=Props('H','P',P(i)/1000,'Q',0,name)*1000;
end
%%

x=1;
P=[P0:5000:Pc];
for i=1:length(P)
    results= main(cellArray,[choix_param, P(i),x]);
    h_1(i)=results(5);
    s_1(i)=results(6);
    T_1(i)=results(2);
    %h_1CP(i)=Props('H','P',P(i)/1000,'Q',1,name)*1000;
end

%%

choix_param= 3; % 1= (v,T); 2= (v,P); 3= (P,T); 4=(T,s); 5= (P,s); 6= (P,h), 7= Ajouter un fluide 

for i=1:length(P)/2
results= main(cellArray,[choix_param, P(i), T_test]);
    h_c(i)=results(5);
    h_cCP(i)=Props('H','P',P(i)/1000,'T',T_test,name)*1000;
end

%%
choix_param= 3; % 1= (v,T); 2= (v,P); 3= (P,T); 4=(T,s); 5= (P,s); 6= (P,h), 7= Ajouter un fluide 

T=[Ttp:1:Tc];
for i=1:length(T)/2
results= main(cellArray,[choix_param, P_test,T(i)]);
    s_c(i)=results(6);
    %s_cCP(i)=Props('S','P',P_test/1000,'T',T(i),name)*1000;
end
%%
figure
plot(h_0/1000,P/1000)
hold on
plot(h_1/1000,P/1000)
 title('P h  ',name)
plot(h_c/1000,P(1:length(h_c))/1000)
plot((h_cCP-h_cCP(1)+h_c(1))/1000,P(1:length(h_cCP))/1000)
%plot(T-273.15,Psat_CP/1000)
ylabel('Pressure [kPa]')
xlabel('Enthalpy [kJ/kg]')
yyaxis right
%plot(h_c/1000,h_c./(h_cCP-h_cCP(1)+h_c(1))*100-100)
%ylabel('Error [%]')
legend('sat curve liquid','sat curve vapor','ADVENS','Coolprop')

figure
plot(s_0,T_0-273.15)
hold on
plot(s_1,T_1-273.15)

 plot(s_c,T(1:length(s_c))-273.15)
 plot(s_cCP-s_cCP(end)+s_c(end),T(1:length(s_cCP))-273.15)
%plot(T-273.15,Psat_CP/1000)
 title('T s  ',name)
ylabel('Temperature [°C]')
xlabel('Entropy [J/K]')
legend('sat curve liquid','sat curve vapor','ADVENS','Coolprop')
%yyaxis right
%plot(s_c,s_c./(s_cCP-s_cCP(end)+s_c(end))*100-100)
%ylabel('Error [%]')

