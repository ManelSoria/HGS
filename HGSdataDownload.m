function HGSdataDownload
%**************************************************************************
%
% HGSdataDownload
%
%**************************************************************************
%
% HGSdataDownload obtains data from BURCAT page and converts it in HGSdata
% database that will be used in the HGS 2.X functions.
% http://garfield.chem.elte.hu/Burcat/THERM.DAT
% HGSdata structure is formed by:
%               .name Contains all the species names
%               .state Contains the state of the species
%               .lim Contains de limits of the Burcat coeficients
%                       lim(ii,1) [K] Lowest temperature limit
%                       lim(ii,2) [K] Temperature where it changes between
%                                       LV and HV
%                       lim(ii,3) [K] Highest temperature limit
%               .spec Contains the elements of this species (Max. 4)
%               .nspec Contains the elemnt number of this species
%               .HV Contains the Burcat coeficients for the temperature
%                       range beetween lim(ii,2) and lim(ii,3)
%               .LV Contains the Burcat coeficients for the temperature
%                       range beetween lim(ii,1) and lim(ii,2)
%               .Mm Contains the molar mass of the species
%               .ena Contains the elements of this species 
%               .nat Contains the elemnt number of this species
%               .comb Contains the names of the mixtures
%               .cspec Contains the species that forms the mixture
%               .cper Contains the species percentage that forms the
%                        mixture
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
%
% Outputs:
%--------------------------------------------------------------------------
% Creates HGSdata.mat
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

clear all;
clc;

global HGSdata

% Visualize name modifications
info = 1;

% Load HGSdata to be updated or create the structure
try
    HGSload
catch
    HGSdata = struct();
end

% This was hard to do
% Do NOT try this at home

download = 1;
wname = 'http://garfield.chem.elte.hu/Burcat/THERM.DAT';
fname = 'DATA_7_coef.txt';

if download
    try
        websave(fname,wname);
    catch 
        error('Uhhh cannot find %s',wname);
    end
end

[fi]=fopen(fname,'r');

if fi==-1
    error('Cannot open input file <%s>',fname);
end

if info 
    fprintf('%s downloaded and saved as %s .. processing \n',wname,fname);
end

% Clean all info data until data start
while 1
    l = fgetl(fi);
    if strcmpi('THERMO ALL',strtrim(l))
       l = fgetl(fi);
       break;
    end
end

ii = 1;
% Loop to obtain data
while 1
    if (ii==1 || mod(ii,100)==0) && info 
        fprintf('Reading ii=%d \n',ii);
    end
    
    l = fgetl(fi);
    
    % End of the data
    if -1 == l
       break;
    end
    
   name = l(1:18);

   if strcmpi(strtrim(name),'Air')==1 % we ignore air as composition is not given
       l = fgetl(fi); l = fgetl(fi); l = fgetl(fi); % skip 3 lines
       continue
   end
   
   if strcmpi(name(1:5),'MgCL2')==1 % we ignore air as composition is not given
       l = fgetl(fi); l = fgetl(fi); l = fgetl(fi); % skip 3 lines
       continue
   end
   
   % Process first line
   % Get name
   HGSdata.name{ii,1} =  name; 
   if strcmp(l(45),'.') % Process exception 
       HGSdata.state{ii,1} =  l(46); % State Triplet
       lcomp = 1;
   else % Usually they begin at 45
       HGSdata.state{ii,1} =  l(45); % State
       lcomp = 0;
   end

   
   % Process temperature limits
   HGSdata.lim(ii,1) = str2double(strtrim(l(49:55)));
   HGSdata.lim(ii,2) = str2double(strtrim(l(67:73)));
   HGSdata.lim(ii,3) = str2double(strtrim(l(58:65)));
   
  
   
   % Process composition
   nn = strfind(l(25+lcomp:44+lcomp),'.');


   for jj=1:4 % There are, at max, 4 elements
       if length(nn) < 4 
            pos = (jj-1)*5 + 25;
            if ~strcmp(l(pos),' ') && ~isstrprop(l(pos),'digit')
                HGSdata.spec{ii,jj} = strtrim(l(pos:pos+1));
                HGSdata.nspec{ii,jj} = str2double(strtrim(l(pos+3:pos+4)));
            end
       else % Normal case
           pos = nn(jj) + 25 + lcomp - 1;
           if jj == 1
               pos2 = 25 + lcomp;
           else
               pos2 = 25 + lcomp + nn(jj-1);
           end

           if ~strcmp(l(pos2),' ')
                HGSdata.spec{ii,jj} = strtrim(l(pos2:pos2+1));
                HGSdata.nspec{ii,jj} = str2double(strtrim(l(pos-2:pos-1)));
           end
       end
   end

   % Rest of the lines, we get the data
   s = zeros(15,1);
   for jj=1:3 % each line
       l = fgetl(fi);
       for q=1:5
           num = l((q-1)*15+1:(q-1)*15+15);
           if ~(q == 5 && jj ==3)
               s((jj-1)*5+q,1) = sscanf(num,'%e');
           end
       end
   end
   HGSdata.HV(ii,:) = s(1:7)';
   HGSdata.LV(ii,:) = s(8:14)';

   ii = ii + 1;


end


%% Cleaning and Reajusting
if info
    fprintf('Changing uppercase in the element names...\n');
end

% Species composition 2nd letter of the elements to lowercase 
% (e.g. AL -> Al )
for ii=1:size(HGSdata.spec,1)
    for jj = 1:size(HGSdata.spec,2)
        if length(HGSdata.spec{ii,jj}) > 1
            HGSdata.spec{ii,jj} = [HGSdata.spec{ii,jj}(1) lower(HGSdata.spec{ii,jj}(2))];
        end
    end
end


if info
    fprintf('Name modifications... ')
end


for ii=1:size(HGSdata.name,1)
    if info
        fprintf('\n <%d> %s',ii,HGSdata.name{ii,1})
    end
    
    % Change name  NE PO BR PB AL and CL 

    % Aluminium AL -> Al
    LowcaseSp('AL')
    
    % Chlorine CL -> Cl
    LowcaseSp('CL')
    
    % Polonium PO -> Po
    LowcaseSp('PO')
    
    % Lead PB -> Pb
    LowcaseSp('PB')
    
    % Bromine BR -> Br
    LowcaseSp('BR')
    
    % Argon AR -> Ar
    LowcaseSp('AR')
    
    % Iridium IR -> Ir
    LowcaseSp('IR')
    
    % Neon NE -> Ne
    LowcaseSp('NE')
    
    % Osmium OS -> Os
    LowcaseSp('OS')
        
    % E- 
    e = strfind(HGSdata.name{ii,1},'electron');
    if ~isempty(e)
       HGSdata.name{ii,1} = strtrim(HGSdata.name{ii,1}(1:e(1)-1));
       % if info; fprintf(' %s -> %s \n',oname,HGSdata.name{ii,1}); end   
    end
    
    
    % Liquid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Liquid change (L) -> (l)
    StateChange('(L)','(l)')
    
    % Liquid change (liq) -> (l)
    StateChange('(liq)','(l)')
    
    % Solid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Solid change (S) -> (s)
    StateChange('(S)','(s)')
    
    % Solid change cr -> (cr)
    StateChange(' cr','(cr)')
    StateChange('(cr)A','(cr.A)')
    StateChange('(cr)B','(cr.B)')
    StateChange('(cr)C','(cr.C)')
    StateChange('(cr)II','(cr.II)')
    StateChange('(cr)I','(cr.I)')
    StateChange('(III)cr','(cr.III)')
    
    
    %  GAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Gas change (G) -> (g)
    StateChange('(G)','(g)')
      
    % Gas change (gas) -> (g)
    StateChange('(gas)','(g)')
    
    % Gas change gas -> (g)
    StateChange(' gas','(g)')
    
    
    % Gas change g -> (g)
    StateChange(' g ','(g)')
    
    % Separte names%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Form = strfind(HGSdata.name{ii,1},'Form');
    if ~isempty(Form)
        if strcmp(HGSdata.name{ii,1}(Form(1)-1),'-')
            HGSdata.name{ii,1} = [HGSdata.name{ii,1}(1:Form(1)-1) ' ' HGSdata.name{ii,1}(Form(1):end)];
            printChange
        end
    end
    
    % Comas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mirar como hacer eso
    HGSdata.name{ii,1} = strtrim(HGSdata.name{ii,1});
%     Comas = strfind(HGSdata.name{ii,1},',');
    blank = strfind(HGSdata.name{ii,1},' ');
%     if ~isempty(Comas) && ~isempty(blank)
%         if Comas(1)<blank(1) 
%             HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},',',' ');
%             printChange
%         end
%      disp('x')  
%     end
           
    %  BLANK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    blank = strfind(HGSdata.name{ii,1},' ');
    if ~isempty(blank)
        
        HGSdata.name{ii,1} =[HGSdata.name{ii,1}(1:blank(1)-1) '[' strtrim(HGSdata.name{ii,1}(blank(1):end)) ']'];
         
        
        % REF ELEMT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eliminate REF ELEMENT and Variants
        StateChange('REF  ELEMENT','')
        StateChange('REF ELEMENT','')
        StateChange('REF  ELEMEN','')
        StateChange('REF ELEMEN','')
        StateChange('REF  ELEME','')
        StateChange('REF ELEME','')
        StateChange('REF  ELEM','')
        StateChange('REF ELEM','')
        StateChange('REF  ELE','')
        StateChange('REF ELE','')
        StateChange('REF  EL','')
        StateChange('REF EL','')
        StateChange('REFELEMENT','')
        StateChange('REF[ELEMENT','[')
        StateChange('REF[ELEME','[')
        StateChange('REFELEM','')
        StateChange('REFEL','')
        StateChange('REFER ELE','')
        StateChange('Ref Elemen','')
        StateChange('REF','')

        % HF298
        StateChange('HF298','')
        
        % JP10
        StateChange('JP10','') 

        %  Radicals Cations and Anions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cat variant
        Cat = strfind(lower(HGSdata.name{ii,1}),'cat');
        if ~isempty(Cat)
            if ~strcmp(HGSdata.name{ii,1}(Cat(1)+3),'e')
                HGSdata.name{ii,1} = [strtrim(HGSdata.name{ii,1}(1:Cat(1)-1)) ']'];
            end
            printChange
        end
        
        % Ani variant
        StateChange('tripleanion','')
        StateChange('bianion','')
        StateChange('anion  H','')
        StateChange('anion','')
        StateChange('Anion','')
        StateChange('anio','')
        StateChange('Anio','')
        
        % Radi variant

        StateChange('biradical','');
        StateChange('Radical','');StateChange('RADICAL','');StateChange('radical','');
        StateChange('Radica','');StateChange('RADICA','');StateChange('radica','');
        StateChange('Radic','');StateChange('RADIC','');StateChange('radic','');
        StateChange('Radi','');StateChange('RADI','');StateChange('Radi','');
        StateChange('Rad','');StateChange('RAD','');StateChange('rad','');
        StateChange('[R]','');
        
        % RRHO
        StateChange('RRHO','')
        StateChange('rrho','')
        
        % ANHARMONIC
        ANHAR = strfind(lower(HGSdata.name{ii,1}),'anhar');
        if ~isempty(ANHAR)
            HGSdata.name{ii,1} = [strtrim(HGSdata.name{ii,1}(1:ANHAR(1)-2)) '(I)'];       
            printChange
        end

        % Atom
        StateChange('Atom','')
        StateChange('atom','')
        StateChange('ATOM','')
        
%%%%%%%%%%%%%%%%%%%% Liquid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Liquid change LIQUID -> (l)
        StateChange('liquid','(l)')
        StateChange('LIQUID','(l)')
        StateChange('Liquid','(l)')
        
        % Liquid change LIQ -> (l)
        StateChange('liq','(l)')
        StateChange('LIQ','(l)')
        StateChange('Liq','(l)')
        
        % Liquid change (l)in a bracket to (l) out of bracket
        OutBracket('(l)')
        
        %  Gas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Gas change GAS -> (g)
        StateChange('GAS','(g)')   
        
        % Gas change (g)in a bracket to (g) out of bracket
        OutBracket('(g)')
        
            
        
        
%%%%%%%%%%%%%%%%%%%% Solid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solid change Solid-* in a bracket to (s.*)out of bracket
        StateChange('(s)[Alpha]','(s.A)')
        StateChange('(s)[Beta]','(s.B)')
        StateChange('(s)[Gama]','(s.C)')
        StateChange('Solid-A','(s.A)')
        StateChange('Solid-B','(s.B)')
        StateChange('Solid-C','(s.C)')
        
        % Solid change SOLID -> (s)
        StateChange('Solid','(s)')
        StateChange('SOLID','(s)')
        StateChange('solid','(s)')
        
        % Solid change (s) (s.A) & co in a bracket to -- out of bracket
        SolChange('(s)','(cr)')
        SolChange('(s)','(s.A)')
        SolChange('(s)','(s.B)')
        SolChange('(s)','(s.C)')
        OutBracket('(cr)')
        OutBracket('(s.A)')
        OutBracket('(s.B)')
        OutBracket('(s.C)')
        OutBracket('(s)')

%%%%%%%%%%%%%%%%%%%%% Names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        
        StateChange('Hydrogen Br','')
        StateChange('G3B3','')
        StateChange('calc','')
        StateChange('<GR>','(s)')
        StateChange('amond','')
        StateChange('Bromotric','')
        StateChange('BrIsocyanogen','')
        StateChange('Bromo','')
        StateChange('Brom','')
        StateChange('11B1','')
        StateChange('FREON 1301','')
        StateChange('11B3','')
        StateChange('MethylCl','')
        StateChange('FC-13','')
        StateChange('FREON-12','')
        StateChange('Formyl Ch','')
        StateChange('CarbonylFlu','')
        StateChange('fluoroIo','')
        StateChange('Tri','')
        StateChange('HBFC-22B1','')
        StateChange('FC-20B1','')
        StateChange('FC-23','')
        StateChange('FC-143A ','')
        StateChange('HCFC-22','')
        StateChange('FC-21','')
        StateChange('Chloroform','')
        StateChange('Isocyanic','')
        StateChange('Cyanic','')
        StateChange('Fulminic','')
        StateChange('Acid','')
        StateChange('Aci','')
        StateChange('Mono','')
        StateChange('Carbon','')
        StateChange('FC-1114','')
        StateChange('PentaFluo','')
        StateChange('BENZENE','')
        StateChange('Toluene','')
        StateChange('Calomel','')
        StateChange('Chloronium','')
        StateChange('Clnitrat','')
        StateChange('Chlo','')
        StateChange('Azidic Acid-d','')
        StateChange('Water-DT','')
        StateChange('NitrAmide','')
        StateChange('Hypoflorous','')
        StateChange('Fluoronium','')
        StateChange('fluorodis','')
        StateChange('Fluor','')
        StateChange('HYDROXYL','')
        StateChange('A 2Sigma+','(I)')
        StateChange('Water-T1','')
        StateChange('perox','')
        StateChange('mercur','')
        StateChange('Iodine','')
        StateChange('Iodid','')
        StateChange('Iridium','')
        StateChange('MgDiboride','')
        StateChange('bromoim','')
        StateChange('Peroxynitr','')
        StateChange('Nitroxyl','')
        StateChange('Amonia','')
        StateChange('Nitro','')
        StateChange('Amo','')
        StateChange('Isodiazene','')
        StateChange('Isodiazene+','')
        StateChange('Isodiazene-','')
        StateChange('HYDRAZINE','')
        StateChange('Hydrazine','')
        StateChange('Hydrazin','')
        StateChange('Azidic ac','')
        StateChange('Sodium','')
        StateChange('Hydro','')
        StateChange('Osmium','')
        StateChange('-D3','')
        StateChange('Phosphine','')
        StateChange('Phosphonium','')
        StateChange('Pmonoxide','')
        StateChange('Phosphor','')
        StateChange('GALENA','')
        StateChange('Paladium','')
        StateChange('Hy','')
        StateChange('Polonium','')
        StateChange('Di','')
        StateChange('amond','Diamond')
        StateChange('Platinum','')
        StateChange('Sulfur','')
        StateChange('Sul','')
        StateChange('diox','')
        StateChange('Sulfate','')
        StateChange('Sulfur3','')
        StateChange('Tritium','')
        StateChange('Technetium','')
        StateChange('Telurium','')
        StateChange('Wurzite','')
        
%%%%%%%%%%%%%%%%%%%%% Sing and Trip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StateChange('singlet','Sing')
        StateChange('SINGLET','Sing')
        StateChange('single','Sing')
        StateChange('tripet','Trip')
        StateChange('triplet','Trip')
        StateChange('triple','Trip')
        StateChange('TRIPLET','Trip')
        StateChange('linear','lin')
        StateChange('doublet','doub')
        StateChange('quartet','quar')
        StateChange('Cyclo','Cy')
        StateChange('Cycl','Cy')
        StateChange('cy','Cy')
        
%%%%%%%%%%%%%%%%%%%%% JANAF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StateChange('JANAF','')
        StateChange('JANAF 65','')
        
 %%%%%%%%%%%%%%%%%%%%% ATcT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StateChange('ATcT C','')
        StateChange('ATcT','')
        
%%%%%%%%%%%%%%%% Eliminate [] empty brackets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StateChange('[R]','');
        StateChange('[  ','[')
        StateChange('[ ','[')
        StateChange('  ]',']')
        StateChange(' ]',']')
        StateChange('[]','')
    end
    
        
        
%%%%%%%%%%%%%%%%%%% REF ELEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % REF-ELEM variant
    StateChange('REF-ELEMENT','')

%     if ~strcmp(oname,HGSdata.name{ii,1})
%         fprintf('<%s> becomes <%s> \n',oname,HGSdata.name{ii,1});
%     end
end


%% Molar Mass calculation
if info
    fprintf(' \n Molar Mass calculation ... \n');
end

Mendeleiev = HGSmendeliev; % Molecular mass of each element

for ii = 1:size(HGSdata.name,1)
    Mm = 0;
    for jj= 1:4
        if ~isempty(HGSdata.spec{ii,jj}) && ~strcmp('E',HGSdata.spec{ii,jj})
            Mm = Mm + HGSdata.nspec{ii,jj} * Mendeleiev.(HGSdata.spec{ii,jj});  
        end
    end
    HGSdata.Mm(ii,1) = Mm;
end

%% Change the format before saving

li=HGSdata.spec; % old format
lin=HGSdata.nspec;
ns=size(HGSdata.Mm,1); 
HGSdata.ena={}; % element names for each of the species
HGSdata.nat={}; % number of atoms of each element, for each of the species
for i=1:ns
    el={};
    nel=[];
    for j=1:4
        if isempty(li{i,j})
            break;
        end
        el{end+1}=li{i,j}; % it may be slow, but if you want it fast don't use Matlab, sinner
        % so please ignore the ** warning
        nel(end+1)=lin{i,j};
    end
    HGSdata.ena{end+1}=el;
    HGSdata.nat{end+1}=nel;
end

% Manel hated my struct so I had to change it
HGSdata = rmfield(HGSdata,'spec');
HGSdata = rmfield(HGSdata,'nspec');

%% Brackets out

roman = {'I' 'II' 'III' 'IV' 'V' 'VI' 'VII' 'VIII' 'IX' 'X' 'XI' 'XII' ...
    'XIII' 'XIV' 'XV' 'XVI' 'XVII' 'XVIII' 'XIX' 'XX'};
HGSdata.nameback = HGSdata.name;
rep = HGSdata.name;
for ii = 1:length(HGSdata.name)
    bracket = strfind(HGSdata.name{ii,1},'[');
    if ~isempty(bracket)
        newname = HGSdata.name{ii,1}(1:bracket(1)-1);
        rep{ii} = newname;
        search = ismember(rep(:),newname);
        if ~all( search == 0)
            num = sum(search);
            if num~= 1
                HGSdata.name(ii,1) = strcat(newname,'(',roman(num),')');
            else
                HGSdata.name{ii,1} = newname; 
            end
        end
    end
end


%% Saving data
    

save('HGSdata.mat','HGSdata');
if info
    fprintf('Saving HGSdata ... \n')
end
%% Nested functions

    % Print changes
    function printChange
        if info
            fprintf(' --> %s',HGSdata.name{ii,1})
            % fprintf(' %s -> %s \n',oname,HGSdata.name{ii,1}) 
        end
    end
    
    function LowcaseSp(sp)
        spc = strfind(HGSdata.name{ii,1},sp);
        if ~isempty(spc)
            Change = [sp(1) lower(sp(2))];
            for q=1:length(spc)
               for jj=1:4
                   if strcmp(Change,HGSdata.spec{ii,jj})
                        HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},sp,Change);
                        printChange 
                        break;
                   end
               end
            end  
        end
    end

    function StateChange(old,new)
        State = strfind(HGSdata.name{ii,1},old);
        if ~isempty(State)
            HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},old,new);
            printChange  
        end
    end

    function OutBracket(state)
        corch = strfind(HGSdata.name{ii,1},state);
        if ~isempty(corch) 
            if length(corch) > 1
                if corch(1) < blank
                    for jj=length(corch):-1:2
                        HGSdata.name{ii,1} = [strtrim(HGSdata.name{ii,1}(1:corch(jj)-1)) strtrim(HGSdata.name{ii,1}(corch(jj)+length(state):end))];
                    end
                else
                    for jj=length(corch):-1:1
                        HGSdata.name{ii,1} = [strtrim(HGSdata.name{ii,1}(1:corch(jj)-1)) strtrim(HGSdata.name{ii,1}(corch(jj)+length(state):end))];
                    end
                    HGSdata.name{ii,1} = [HGSdata.name{ii,1}(1:blank(1)-1) state HGSdata.name{ii,1}(blank(1):end)];
                end
                printChange  
            elseif corch > blank(1) 
                name = [HGSdata.name{ii,1}(1:blank(1)-1) state];
                comp = [strtrim(HGSdata.name{ii,1}(blank(1):corch(1)-1)) strtrim(HGSdata.name{ii,1}(corch(1)+length(state):end))];
                HGSdata.name{ii,1} = [name comp];
                printChange  
            end
        end
    end

    function SolChange(old,new)
        corch1 = strfind(HGSdata.name{ii,1},old);
        corch2 = strfind(HGSdata.name{ii,1},new);
        if ~isempty(corch1) &&  ~isempty(corch2)
            if corch1(1)<corch2(1)
                HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},new,'');
                HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},old,new);
                printChange
            else
                HGSdata.name{ii,1} = strrep(HGSdata.name{ii,1},old,'');
                printChange
            end
        end
    end

end

