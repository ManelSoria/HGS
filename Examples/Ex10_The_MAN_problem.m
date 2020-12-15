%***********************************************************************************************************
%* HGS 2.0
%* Original Caleb Fuster
%
%* ESEIAAT UPC          
%***********************************************************************************************************
%
% Example 10: The MAN problem
% Wet air reaction with hydrocarbons
clear all
clc

hum = 0:2:100;

for ii = 1:length(hum)
    % Humidity ratio
    x = hum(ii); % [%] x = Kgwater/Kgdryair *100

    % Wet air in Kg
    Kgwetair = 0.1; % [kg] Kgwetair = Kgwater + Kgdryair
    Kgdryair = Kgwetair / (x/100 + 1); % [kg]
    Kgwater = x / 100 * Kgdryair; % [kg]

    % Dry air composition in mols
    dryaircomp = [0.79 0.21];

    % Molar mass
    Mmdryair = HGSprop({'N2' 'O2'},dryaircomp,[],[],'Mm');
    Mmwater = HGSprop({'H2O'},1,[],[],'Mm');

    % Mols air
    nwater = Kgwater * 1000 / Mmwater;
    ndryair = Kgdryair * 1000 / Mmdryair;

    nN2dry = dryaircomp(1) * ndryair;
    nO2dry = dryaircomp(2) * ndryair;

    nwetaircomp = [nN2dry nO2dry nwater] / (nN2dry + nO2dry + nwater);

    % Wet air
    wet = {'N2' 'O2' 'H2O'};
    nwetair = Kgwetair * 1000  / HGSprop(wet,nwetaircomp,[],[],'Mm'); % [mols] 
    nwet =  nwetair * nwetaircomp; % [mols]

    % Fuels
    fuel = {'CH4'};
    nfuel = [0.4];  % [mols]

    % Species
    % HGS functions requires all the species that you want to be calculated.
    % If you want to put dissociation, it is possible. Like in this example H2
    % and H or O2 and O.

    species = {wet{1:end} fuel{1:end} 'H2' 'CO2' 'NO2'};
    n0  =     [   nwet        nfuel     0    0     0  ];

    % Suposing that Wet air and Fuels have not the same inlet T
    Hwet = HGSprop(wet,nwet,298.15,1,'H');
    Hfuels = HGSprop(fuel,nfuel,350,1,'H');

    Hin = Hwet + Hfuels;

    % HGStp

    [Tp(ii),n(ii,:),species,~] = HGStp(species,n0,'H',Hin,1);
    nper(ii,:) = n(ii,:)/sum(n(ii,:));
    
end

figure
plot(hum,Tp);

figure
hold on
for ii = 1 : length(species)
    plot(hum,nper(:,ii));
end
legend(species)
