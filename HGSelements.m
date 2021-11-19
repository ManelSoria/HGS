function [eln,elq] = HGSelements(name)
%**************************************************************************
%
% HGSelements(name)
%
%**************************************************************************
%
% HGSelements returns the elements present in the species name and their
% number
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name or code of the species
%
%**************************************************************************
% Outputs:
%--------------------------------------------------------------------------
% eln --> List with the name of the elements present in the species
% elq --> Vector with their quantities
%
%**************************************************************************
% Examples:
% [eln,elq]=HGSelements('H2O');
% [eln,elq]=HGSelements(2247);
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global HGSdata;HGSload;

id=HGSid(name);

if length(id)>1
   error('HGSelements works 1 by 1 species. Do not enter more than 1 species') 
end

eln=HGSdata.ena{id};
elq=HGSdata.nat{id};
end