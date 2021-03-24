function [Res] = HGSsingle(species,property,T,P)
%**************************************************************************
%
% [Res] = HGSsingle(species,property,T,P)
%
%**************************************************************************
%
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% property --> Property requested (see below)
% T --> [K] Temperature
% P --> [bar] Pressure
%
% Outputs:
%--------------------------------------------------------------------------
% Res --> Property result
%           Mm [g/mol]
%           Cp [kJ/(mol*K)]
%           Cv [kJ/(mol*K)]
%           h [kJ/mol]
%           s [kJ/(mol*K)]
%           g [kJ/mol]
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    


code = zeros(7,1); 
% code: (num)- Property (need to be calculated before)
%        (1) - Molar mass  ()
%        (2) - cp          (7- Burcat)
%        (3) - cv          (7- Burcat & 2- Cp)
%        (4) - h           (7- Burcats)
%        (5) - s           (7- Burcats)
%        (6) - g           (7- Burcats & 4- H & 5- S)
%        (7) - Burcat Coef ()


switch property
    case 'Mm'
        code(1) = 1;
    case 'cp'
        code(2) = 1;
        code(7) = 1;
    case 'cv'
        code(2:3) = 1;
        code(7) = 1;
    case 'h'
        code(4) = 1;
        code(7) = 1;
    case 's'
        code(5) = 1;
        code(7) = 1;
    case 'g'
        code(4:7) = 1;
    otherwise 
        error('Uknown property %s ',property);
end

id = HGSid(species);
global HGSdata;
HGSload;

%% Burcat coeficients
if code(7)  
    if T <= HGSdata.lim(id,1) || T >= HGSdata.lim(id,3)
        lim(1,1) = HGSdata.lim(id,1);
        lim(1,3) = HGSdata.lim(id,3);
        name = species;
        error('HGSsingle: Ups... Temperature is not between the limits (%.2fK-%.2fK) for %s',lim(1,1),lim(1,3),name)
    end  

    if T <= HGSdata.lim(id,2)
        a(1,:) = HGSdata.LV(id,:);
    else
        a(1,:) = HGSdata.HV(id,:);
    end 
end

%% Molar Mass
if code(1)
     Res = HGSdata.Mm(id);
end

%% Cp
if code(2)
    [cp] = HGScp(a,T);
    
    Res = cp;  
end

%% Cv
if code(3)
    [Res] = HGScv(cp);
end

%% H
if code(4)
    [h] = HGSh(a,T); 

    Res = h;
end

%% S
if code(5)
    [s] = HGSs(a,T,P,HGSdata.state{id});
    Res = s;
end

%% G
if code(6)
    [Res] = HGSg(s,h,T);
end


end