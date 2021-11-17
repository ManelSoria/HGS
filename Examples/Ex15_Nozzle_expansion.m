%***********************************************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC      
%***********************************************************************************************************
%
% Nozzle expansion process.
% LH2-LOX reaction

% Exemplification of the differents variables across of a nozzle expansion
% Changing the for of the P vector by a solver is what HGSnozzle does.

close all
clear

% The inlet conditions to combustion chamber of the Vinci engine are:
species={'H2','O2','H2O','OH','O','H'};
n0 = [2875.496,1052.566,0,0,0,0]; % mol/s (obtained from kg/s)
P0=62;
Pa=0; % Vacuum
Hin=-19.0609; % kJ (liquid inlet, obtained with INIST)

[Tp,~,np,~] = HGStp(species,n0,'H',Hin,P0);

%% Nozzle expansion
% We generate a vector with the pressure points to be obtained
% From 10 bar below the chamber pressure to 0.01 in 2 diferent steps
% 10 bar is arbitrary, we just want to make sure that we begin before the
% throat but far from the inlet, where the velocity would be 0 and the area
% infinite.
P = [P0-10:-2:1, linspace(1,0.01,20)]; 
[species,n,T,v,M,A,F,Isp] = HGSnozzle(species,np,Tp,P0,P,Pa,'Shifting');

% Plots

figure
plot(P,T, '-b','LineWidth',2)
set ( gca, 'xdir', 'reverse')
title('T versus P')
xlabel('P [bar]')
ylabel('T [K]')
set(gca,'FontSize',18)

%%
figure
hold on
plot(P,n(3,:)./sum(n(:,3)),'LineWidth',2)
set ( gca, 'xdir', 'reverse' )
title('H2O fraction versus P')
xlabel('P [bar]')
ylabel('H2O Molar fraction []')
set(gca,'FontSize',18)
%%
figure
plot(P,v,'LineWidth',2)
title('v vs P')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('v [m/s]')
set(gca,'FontSize',18)
%%
figure
plot(P,M,'LineWidth',2)
hold on
plot(P,ones(size(P)),'r')
title('Mach vs P')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('Mach []')
set(gca,'FontSize',18)

%%
% We can find the exact throat pressure with:
[~,~,~,Pt,flag] = HGSisentropic(species,np,Tp,P0,'Shifting','M',1);
plot([Pt Pt],[0 6],'r')


%%
figure
plot(P,F,'LineWidth',2)
title('Thrust vs P')
set ( gca, 'xdir', 'reverse' )
xlabel('P [bar]')
ylabel('F [N]')
set(gca,'FontSize',18)

%%
figure
D=1000*2*sqrt(A/pi);
semilogy(P,D,'LineWidth',2)
hold on
semilogy([Pt Pt],[min(D)/2 max(D)],'r')
set ( gca, 'xdir', 'reverse' )
title('Diameter vs P')
xlabel('P [bar]')
ylabel('Diameter [mm]')
set(gca,'FontSize',18)

Dt=interp1(P0-P,D,P0-Pt); % P0-P to have x values in ascending order
De=D(end);

fprintf('Throat diameter %.2f mm\nExit diameter %.2f mm \n',Dt,De);
fprintf('Nozzle expansion ratio = %.2f \n',(De/Dt)^2);
fprintf('Thrust = %.2f kN\n',F(end)/1000);