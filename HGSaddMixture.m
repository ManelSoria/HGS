function  HGSaddMixture(name,species,percent)
%**************************************************************************
%
% HGSaddMixture(name,species,percent)
%
%**************************************************************************
%
% HGSaddMixture add a mixture to the HGSdata database
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name of the mixture f.e. 'Air'
% species --> species from the mixture f.e. {,...}
% percent --> Percentage of each species in the mixture [,...]
%
% Outputs:
%--------------------------------------------------------------------------
% HGSdata update
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global HGSdata; HGSload

% Check if no species or mixture use the name
try
    [~] = HGSid(name);  
    flag = 1;
catch
    flag = 0;
end

if flag
   error('Ups,.. this name is already used') 
end

% Check if all species are in HGSdata database
try
    [~] = HGSid(species);
catch
   error('Ups,.. at least one of the species is not in the HGSdata')
end

if isnumeric(species)
    for ii=1:length(species)
        spec{ii} = HGSdata.name{species(ii)};
    end
    species = spec;
end

if ~(sum(percent) == 100)
   error('Ups,... The total percentage is not 100%%') 
end

if ~ischar(name)
    name = name{1};
end

% If the HGSdata has no mixtures
if ~isfield(HGSdata,'comb')
    HGSdata.comb = {};
    HGSdata.cspec = {};
    HGSdata.cper = {};
end

id = HGSid(species);
lm = length(HGSdata.Mm);

if max(id) > lm
    buildS = {};
    buildP = [];
   for ii = 1:length(id)
       if id(ii) > lm
           for jj = 1:length(HGSdata.cspec{id(ii)-lm})
               buildS{end+1} = HGSdata.cspec{id(ii)-lm}{jj};
               buildP(end+1) = HGSdata.cper{id(ii)-lm}(jj)*percent(ii)/100;
           end
           
       else
           buildS{end+1} = species{ii};
           buildP(end+1) = percent(ii);
       end
   end  
   species = buildS;
   percent = buildP;
end

HGSdata.comb{end+1,1} = name;
HGSdata.cspec{end+1,1} = species;
HGSdata.cper{end+1,1} = percent;

end