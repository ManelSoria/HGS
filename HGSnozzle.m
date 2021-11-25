function  [species,n,T,v,M,A,F,Isp] = HGSnozzle(species,n0,T0,P0,P,Pa,Fro_Shift)
%**************************************************************************
%
% [species,n,T,v,M,A,F,Isp] = HGSnozzle(species,n0,T0,P0,Pe,Pa,Fro_Shift,
%                             options1,options2)
%
%**************************************************************************
% 
% HGSnozzle evaluates different flow properties as a function of the
% pressure, during a isentropic expansion beginning with a very low
% velocity
%  
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or code of inlet species
% n0 --> [mols] Number of mols/s of each inlet species
% T0 --> [K] Inlet temperature
% P0 --> [bar] Inlet pressure
% P -->  [bar] Pressure vector (with all the pressures that have to be
%         evaluated)
% Pa --> [bar] Atmospheric pressure
% Fro_Shift --> Select between: 'Frozen' for frozen flow
%                               'Shifting' for shifting flow
%
% Outputs:
%--------------------------------------------------------------------------
% species --> String or code of species
% n --> [mols] Matrix of pecies mols, sorted as: n(species, pressure)
% T --> [K] Exit temperature
% M --> [] Exit Mch
% A --> [m^2] Exit area
% F --> [N] Thrust
% Isp --> [s^]Specific impulse, g0 = 9.807 m/s^2
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

[id] = HGSid(species);

% Total mass
[mm] = HGSprop(id,n0,T0,P0,'Mm'); % g/mol
m=sum(n0)*mm*1e-3; % kg/s
g0 = 9.807;

if ~strcmp(Fro_Shift,'Shifting') && ~strcmp(Fro_Shift,'Frozen')
    error('Your variable Fro_Shift is no one accepted by this function. Only Frozen and Shifting are accepted')
end

for ii=1:length(P)
    fprintf('P = %f,  %i /%i \n',P(ii),ii,length(P))
    [T(ii),~,n(:,ii),M(ii),flag] = HGSisentropic(id,n0,T0,P0,Fro_Shift,'P',P(ii)); % K , mols, [],...
    if flag ~=1, error('HGSnozzle failed to converge/1 flag=%d',flag), end  
    [Rg,a(ii),Mm] = HGSprop(id,n(:,ii),T(ii),P(ii),'Rg','a','Mm'); % kJ/(kg*K), m/s, g/mol
    rho = P(ii)*1e5/(Rg*1000*T(ii)); % kg/m^3 Convert bar to Pa g 2 kG
    v(ii)=M(ii)*a(ii); %m/s
    A(ii) = m/(v(ii)*rho); % m^2
    F(ii) = m*v(ii)+A(ii)*(P(ii)-Pa)*1e5; % Convert bar to Pa
    Isp(ii) = F(ii)/(m*g0);
end


end