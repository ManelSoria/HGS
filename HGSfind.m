function  HGSfind(name,complete)
%**************************************************************************
%
% HGSfind(name,complete)
%
%**************************************************************************
%
% HGSfind finds the species and the mixtures that contain the string name.
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> String to be found
% complete --> Prints the complete name if complete = 1
%
% Outputs:
%--------------------------------------------------------------------------
% Command Window Print
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    


global HGSdata; HGSload
maxid = length(HGSdata.Mm);

if ~exist('complete','var') || complete == 0
    Name_base = HGSdata.name;
elseif   complete == 1
    Name_base = HGSdata.nameback;
else
    error('Ups,... Check complete parameter (Only allowed 0 and 1)')
end


cap = 1; % in case there is no name with the string
fprintf('Species that contain %s \n',name)
for ii=1:maxid
   if strfind(Name_base{ii},name)
       cap = 0;
       fprintf('<%d>  %s \n',ii,Name_base{ii})
   end
end
if cap % Print None to be a fancy list
    fprintf('None \n')
end


cap = 1; % in case there is no name with the string
if isfield(HGSdata,'comb') && ~isempty(HGSdata.comb)
    fprintf('Mixtures that contain %s \n',name)
    for ii=1:length(HGSdata.comb)
       if strfind(HGSdata.comb{ii},name)
           cap = 0;
           fprintf('<%d>  %s \n',maxid+ii,HGSdata.comb{ii})
       end
    end
end
if cap % Print None to be a fancy list
    fprintf('None \n')
end

end