% HGSdataRebuild regenerates all the data abd. This is the function to call if
% you suspect problems with your data base file.

clear all;

[p,n,e]=fileparts(matlab.desktop.editor.getActiveFilename);
if ismac || isunix % decent unix
    delim='/';
else % windows
    delim='\';    
end
databasefile=sprintf('%s%sHGSdata.mat',p,delim);

ok=input(sprintf('Ok to erase and rebuild HGSdata located in <%s> (Y/N) ? \n',databasefile),'s');

if ok~='Y'
    fprintf('Bye\n');
    return
end

if exist(databasefile,'file')
    delete(databasefile);
end

HGSdataDownload(databasefile); % Downloads and saves BURCAT data base

% add mixtures

HGSaddMixture('RP1'    ,{'C10H22', 'C10H18'},[61.25, 38.75]); fprintf('RP1 appended \n');
HGSaddMixture('Air7921',{'N2', 'O2'},[79, 21]); fprintf('Air7921 appended \n');

global HGSdata

save(databasefile,'HGSdata');
fprintf('Saved data to <%s> \n',databasefile)
