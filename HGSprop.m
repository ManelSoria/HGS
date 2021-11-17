function [varargout] = HGSprop(species,n,T,P,varargin)
%**************************************************************************
%
% [varargout] = HGSprop(species,n,T,P,varargin)
%     e.g. [Mm,Cp,Cv,H,S,G,Rg,gamma,a] = HGSprop({'H2' 'O2'},[1 2],300,1)
%          [H,S] = HGSprop({'H2' 'O2'},[1 2],300,1,'H','S')
%          [S,H] = HGSprop({'H2' 'O2'},[1 2],300,1,'S','H')
%
%**************************************************************************
%
% HGSprop returns the properties of a mixture of gasses
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or code of species
% n --> [mols] Number of mols per species
% T --> [K] Temperature
% P --> [bar] Pressure
% varargin --> Expected return: 'Mm' 'Cp' 'Cv' 'H' 'S' 'G' 'Rg' 'gamma' 'a'
%                                If it is empty, all the properties will be
%                                returned
%
% Outputs:
%--------------------------------------------------------------------------
% varargout --> Property result
%                Mm [g/mol]
%                Cp [kJ/K]
%                Cv [kJ/K]
%                H [kJ]
%                S [kJ/K]
%                G [kJ]
%                Rg [kJ/(kg*K)]
%                gamma 
%                a [m/s]
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

%% Checks

Inp = length(varargin);

if length(n) ~= length(species)
     error('Ups..., Species and mols lengths are not the same. Check it')
end


global HGSdata; HGSload
global R; HGSr

[id] = HGSid(species);

% Rebuild mixtures
if max(id) > length(HGSdata.Mm)
   [species,n] = HGSrebuild(species,n);
   [id] = HGSid(species);
end


%% function


code = zeros(10,1); 

% code: (num)- Property (need to be calculated before)
%        (1) - Molar mass  ()
%        (2) - Cp          (10- Burcat)
%        (3) - Cv          (10- Burcat & 2- Cp)
%        (4) - H           (10- Burcats)
%        (5) - S           (10- Burcats)
%        (6) - G           (10- Burcats)
%        (7) - Rg          (1- Mm)
%        (8) - Gamma       (2- Cp & 3- Cv)
%        (9) - a           (1- Mm & 2- Cp & 3- Cv & 7- Rg & 8- Gamma)
%        (10)- Burcat Coef ()

if isempty(varargin)
    Inp = 1;
    code(1:10) = 1;   
else 
    for ii=1:length(varargin)
        switch varargin{ii}
            case 'Mm'
                code(1) = 1;
            case 'Cp'
                code(2) = 1;
                code(10) = 1;
            case 'Cv'
                code(1:3) = 1;
                code(10) = 1;
            case 'H'
                code(4) = 1;
                code(10) = 1;
            case 'S'
                code(5) = 1;
                code(10) = 1;
            case 'G'
                code(4:6) = 1;
                code(10) = 1;
            case 'Rg'
                code(1) = 1;
                code(7) = 1;
            case 'gamma'
                code(2:3) = 1;
                code(8) = 1;
                code(10) = 1;
            case 'a'
                code(1:3) = 1;
                code(7:10) = 1;

            otherwise
                error('Ups,... %s is not a property is not implemented in HGSprop',varargin{ii})
        end
    end

end

spec = length(species);

%% Burcat coeficients

if code(10)    
    if length(T) ~= 1
         error('Ups..., Temperatures length has multiple options. Check it')
    end
    for ii=1:spec
        if T < HGSdata.lim(id(ii),1) || T > HGSdata.lim(id(ii),3)
           lim(1,1) = HGSdata.lim(id(ii),1);
           lim(1,3) = HGSdata.lim(id(ii),3);
           
           name = id(ii);
           error('HGSsingle: Ups... Temperature (%f) is not between the limits (%f - %f) for %d.(Use HGSfind to know the element)',...
                                                  T                        ,lim(1,1),lim(1,3)  ,name)
        end  

        if T <= HGSdata.lim(id(ii),2)
            a(ii,:) = HGSdata.LV(id(ii),:);
        else
            a(ii,:) = HGSdata.HV(id(ii),:);
        end 
    end
end

%% Molar Mass
if code(1)
    Mm_i = zeros(spec,1);

    for ii=1:spec
        Mm_i(ii) = HGSdata.Mm(id(ii));
    end
   
    Mm = dot(n,Mm_i)/sum(n);
end

%% Cp
if code(2)
    cp_i = zeros(spec,1);
     for ii=1:spec
         [cp_i(ii)] = HGScp(a(ii,:),T);
     end
     cp = dot(n,cp_i);
end

%% Cv
if code(3)
    for ii=1:spec
         [cv_i(ii)] = HGScv(cp_i(ii));
     end
     cv = dot(n,cv_i);
end

%% H
if code(4)
    h_i = zeros(spec,1);
    for ii=1:spec
        h_i(ii) = HGSh(a(ii,:),T);
    end
    h = dot(n,h_i);
end

%% Partial pressure for S and G calculations
if code(5)
    nt = 0;
    for ii=1:spec
        if strcmp(HGSdata.state{id(ii)},'G')
            nt = nt + n(ii);
        else
            error('Ups,.. Right now entropy can be calculated only for gas mixtures')
        end
    end
    
    P_i = P.*n/nt;
end

%% S
if code(5)
    s_i = zeros(spec,1);
    
    for ii=1:spec
        s_i(ii) = HGSs(a(ii,:),T,P_i(ii),HGSdata.state{id(ii)});
    end
    
    s = dot(n,s_i);
end

%% G
if code(6)
    g_i = zeros(spec,1);
    for ii=1:spec 
        g_i(ii) = HGSg(s_i(ii),h_i(ii),T);
    end
    g = dot(n,g_i);
end

%% Rg
if code(7)
    Rg = R/Mm*1000;
end

%% Gamma
if code(8)
    gamma = cp/cv;
end

%% a
if code(9)
   sound = sqrt(gamma*Rg*1e3*T);
end

%% Outputs
if isempty(varargin)
    varargout{1} = Mm;
    varargout{2} = cp;
    varargout{3} = cv;
    varargout{4} = h;
    varargout{5} = s;
    varargout{6} = g;
    varargout{7} = Rg;
    varargout{8} = gamma;
    varargout{9} = sound;
else
    for ii=1:Inp
        switch varargin{ii}
            case 'Mm'
                varargout{ii} = Mm;
            case 'Cp'
                varargout{ii} = cp;
            case 'Cv'
                varargout{ii} = cv;
            case 'H'
                varargout{ii} = h;
            case 'S'
                varargout{ii} = s;
            case 'G'
                varargout{ii} = g;
            case 'Rg'
                varargout{ii} = Rg;
            case 'gamma'
                varargout{ii} = gamma;
            case 'a'
                varargout{ii} = sound;

        end   
   end
end



