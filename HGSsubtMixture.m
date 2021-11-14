function  HGSsubtMixture(name)
%**************************************************************************
%
% HGSsubtMixture(name)
%
%**************************************************************************
%
% HGSsubtMixture eliminates a mixture from the HGSdata database
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name of the mixture f.e. 'Air'
% 
% Outputs:
%--------------------------------------------------------------------------
% HGSdata update
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global HGSdata; HGSload


% Check if no species or mixture use the name
try
    [id] = HGSid(name);
catch
    error('Ups,.. this name is not used in HGSdata')
end

del = id - length(HGSdata.Mm);

HGSdata.comb(del) = [];
HGSdata.cspec(del) = [];
HGSdata.cper(del) = [];

end