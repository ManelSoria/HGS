%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC      
%***********************************************************************************************************
%
% The MAN problem
% Wet air reaction with hydrocarbons. Using Relative humidity.
% At ambient pressure and temperature.

clearvars
close all
clc


RH = 0:2:100; % [%] Relative humidity
T1 = 298.15; % [K] Air temperature
T2 = 350; % [K] Fuel temperature
Kgwetair = 1; % [kg] Wet air mass flow
P = 1; % [bar] Atmospheric pressure
V = 0.01; % [m^3] Volume
R = 8.3144621; % [J / (mol*K)]

nwetair = P*10^5*V / (R*T1); % [mol]


% Good approx up to 40ºC. Be carefull Tc (ºC)
% http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/relhum.html#c4
VDsat = @(Tc) 5.018 + 0.32321*Tc + 8.1847*10^(-3)*Tc^2 +3.1243*10^(-4)*Tc^3;
VDsatAmb = VDsat(T1-273.15);


for ii = 1:length(RH)
    % RH = VD / VDsat *100;
    VD = RH(ii) * VDsatAmb / 100 ;
    
    % VD = KgH2O / KgH2 =  KgH2O / 2
    KgH2O = VD * 2;
        
    % Dry air composition in mols
    dryaircomp = [0.79 0.21];
    
    % Molar mass
    Mmdryair = HGSprop({'N2' 'O2'},dryaircomp,[],[],'Mm');
    MmH2O = HGSprop({'H2O'},1,[],[],'Mm');
    
    nH2O = KgH2O / MmH2O;
            
    ndryair  = nwetair - nH2O;
    
    % Wetair composition in mols
    nwetaircomp = [dryaircomp(1)*ndryair dryaircomp(2)*ndryair  nH2O] / (ndryair + nH2O) ;


    % Wet air
    wet = {'N2' 'O2' 'H2O'};
    nwetair = Kgwetair * 1000  / HGSprop(wet,nwetaircomp,[],[],'Mm'); % [mols] 
    nwet =  nwetair * nwetaircomp; % [mols]

    % Fuels
    fuel = {'CH4' 'C2H6' 'C3H8'};
    nfuel = [1.0 0.8 0.5];  % [mols]

    % Species
    % HGS functions requires all the species that you want to be calculated.
    % If you want to put dissociation, it is possible. Like in this example H2
    % and H or O2 and O.

    species = {wet{1:end} fuel{1:end}  'CO2'  'CO'};
    n0  =     [   nwet        nfuel      0     0 ];

    % Suposing that Wet air and Fuels have not the same inlet T
    Hwet = HGSprop(wet,nwet,T1,1,'H');
    Hfuels = HGSprop(fuel,nfuel,T2,1,'H');

    Hin = Hwet + Hfuels;

    % HGStp

    [Tp(ii),species,n(ii,:),~] = HGStp(species,n0,'H',Hin,1);
    nper(ii,:) = n(ii,:)/sum(n(ii,:));
    
end

figure
plot(RH,Tp);
ylabel('Tp [K]')
xlabel('Humidity [%]')

figure
hold on
for ii = 1 : length(species)
    plot(RH,nper(:,ii)*100);
end
xlabel('Humidity [%]')
ylabel('Mols [%]')
axis([0 100 0 25])
legend(species)

figure
hold on
for ii = 1 : length(species)
    plot(RH,n(:,ii));
end
xlabel('Humidity [%]')
ylabel('Mols [mols]')
axis([0 100 0 5])
legend(species)
