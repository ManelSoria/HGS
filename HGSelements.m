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
% Outputs:
% eln --> List with the name of the elements present in the species
% elq --> Vector with their quantities
%--------------------------------------------------------------------------
% Examples:
% [a,b]=HGSelements('H2O');
% [a,b]=HGSelements(2247);
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

if isnumeric(name)
    id=name;
else
    id=HGSid(name);
end

global HGSdata;HGSload;

eln=HGSdata.ena(id);
eln=eln{1};
elq=HGSdata.nat(id);
elq=elq{1};
end