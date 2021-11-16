%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC      
%***********************************************************************************************************
%
% Nozzle expansion process.
% LH2-LOX reaction

% Exemplification of the differents variables across of a nozzle expansion
% Changing the for of the P vector by a solver is what HGSnozzle does.

close all
%% First of all, we need the combustion species :D.
% Look on example 12 for more info about this. Vinci cycle

hH2_INIST=INIST('H2','h_pt',77.5,223.824)*HGSsingle('H2','Mm')/1000; % kJ/mol
hH2_ref_INIST=INIST('H2','h_pt',77.5,350) * HGSsingle('H2','Mm')/1000; % kJ/mol
deltaH_H2=hH2_ref_INIST - hH2_INIST; 
hH2_ref_HGS=HGSsingle('H2','h',350,77.5);
hH2_inlet_HGS = hH2_ref_HGS - deltaH_H2; % kJ/mol

hO2_INIST=INIST('O2','h_pt',77.5,94.518)*HGSsingle('O2','Mm')/1000; % kJ/mol
hO2_ref_INIST=INIST('O2','h_pt',77.5,350) * HGSsingle('O2','Mm')/1000; % kJ/mol
deltaH_O2=hO2_ref_INIST - hO2_INIST; 
hO2_ref_HGS=HGSsingle('H2','h',350,77.5);
hO2_inlet_HGS = hO2_ref_HGS - deltaH_O2; % kJ/mol


species = {'H2', 'O2', 'H2O', 'OH', 'O', 'H'};
n0 = [2875.49603174603,1052.56578536159,0,0,0,0];
Hin = n0(1)*hH2_inlet_HGS + n0(2)*hO2_inlet_HGS;
P0=62;

[Tp,~,np,~] = HGStp(species,n0,'H',Hin,P0);

%% The expansion
% bar 60 to 0.01 in 2 diferent linespace
P = [P0-2:-2:1, linspace(1,0.01,20)]; 
[species,n,T,v,M,A,F,Isp] = HGSnozzle(species,np,Tp,P0,P,Pa,'Shifting');

%% Plots

figure(1)
plot(P,T)
set ( gca, 'xdir', 'reverse' )
title('P vs T')
xlabel('P [bar]')
ylabel('T [K]')

% figure(2)
% hold on
% plot(P,n(1,:)./sum(n(:,1))*100)
% plot(P,n(2,:)./sum(n(:,2))*100)
% plot(P,n(3,:)./sum(n(:,3))*100)
% plot(P,n(4,:)./sum(n(:,4))*100)
% plot(P,n(5,:)./sum(n(:,5))*100)
% plot(P,n(6,:)./sum(n(:,6))*100)
% set ( gca, 'xdir', 'reverse' )
% title('P vs mols')
% xlabel('P [bar]')
% ylabel('% []')
% legend(species)


figure(3)
plot(P,M)
title('Mach vs T')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('Mach []')

figure(4)
plot(P,v)
title('v vs T')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('v [m/s]')

figure(5)
plot(P,F)
title('F vs T')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('F [N]')

figure(6)
plot(P,A)
set ( gca, 'xdir', 'reverse' )
title('P vs A')
xlabel('P [bar]')
ylabel('A [m^2]')

figure(7)
plot(P,sqrt(A/pi))
title('r vs T')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('r [m]')

figure(8)
plot(P,Isp)
title('P vs Isp')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('Isp []')