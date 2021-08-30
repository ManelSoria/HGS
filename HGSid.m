function [id] = HGSid(species)
%**************************************************************************
%
% [id] = HGSid(species)
%
%**************************************************************************
% 
% HGSid finds the id of the species to improve the velocity of the code
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
%
% Outputs:
%--------------------------------------------------------------------------
% id --> Species id
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    


if isnumeric(species)
    id=species;
    return
end

global HGSdata;HGSload;


ns=size(HGSdata.Mm,1); % Number of species
if isfield(HGSdata,'comb') && ~isempty(HGSdata.comb)
    combination = 1;
else
    combination = 0;
end

if iscell(species)    
    id = zeros(1,length(species));

    for ii = 1:length(species)
        id(ii) = find1(species{ii});
    end

    return
end

if ischar(species)
    id = find1(species);
    return
end


function id = find1(name)
    search = ismember(HGSdata.name(:),name);
    if combination
        search1 = ismember(HGSdata.comb(:),name);
    else
        search1 = 0;
    end

    if all(search == 0) && all(search1 == 0)
        error('HGSid:  %s not found in the data base\n',name);
    elseif ~all(search == 0)
        id = find(search);
    else
        id = ns + find(search1);
    end
end
error('uhh ? HGSid wrong data type');



end