function  HGSfind(name,complete)
%**************************************************************************
%
% HGSfind(name,complete)
%
%**************************************************************************
%
% HGSfind finds the species and the mixtures that contains the (name) in
% his name.
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name to find
% complete --> Returns the complete name if complete = 1; No complete = 0.
%
% Outputs:
%--------------------------------------------------------------------------
% Comand Window Print
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    


global HGSdata; HGSload
maxid = length(HGSdata.Mm);

if ~exist('complete','var') || complete == 0
    Names = HGSdata.name;
elseif   complete == 1
    Names = HGSdata.nameback;
else
    error('Ups,... Check complete parameter (Only allowed 0 and 1)')
end

fprintf('Species that contain %s \n',name)
for ii=1:maxid
   if strfind(Names{ii},name)
       fprintf('<%d>  %s \n',ii,Names{ii})
   end
end

if isfield(HGSdata,'comb') && ~isempty(HGSdata.comb)
    fprintf('Mixtures that contain %s \n',name)
    for ii=1:length(HGSdata.comb)
       if strfind(HGSdata.comb{ii},name)
           fprintf('<%d>  %s \n',maxid+ii,HGSdata.comb{ii})
       end
    end
end

end