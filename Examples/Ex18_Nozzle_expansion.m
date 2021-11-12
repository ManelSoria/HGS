%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC      
%***********************************************************************************************************
%
% Nozzle expansion process.
% LH2-LOX reaction

% Exemplification of the differents variables across of a nozzle expansion
% Changing the for of the P vector by a solver is what HGSnozzle does.

%% First of all, we need the combustion species :D.
% Look on example 12 for more info about this.

hH2_INIST=INIST('H2','h_pt',77.5,223.824)*HGSsingle('H2','Mm')/1000; % kJ/mol
hH2_ref_INIST=INIST('H2','h_pt',77.5,350) * HGSsingle('H2','Mm')/1000; % kJ/mol
deltaH_H2=hH2_ref_INIST - hH2_INIST; 
hH2_ref_HGS=HGSsingle('H2','h',350,77.5);
hH2_inlet_HGS = hH2_ref_HGS - deltaH_H2; % kJ/mol

hO2_INIST=INIST('O2','h_pt',77.5,94.518)*HGSsingle('O2','Mm')/1000; % kJ/mol
hO2_ref_INIST=INIST('O2','h_pt',77.5,350) * HGSsingle('O2','Mm')/1000; % kJ/mol
deltaH_O2=hO2_ref_INIST - hO2_INIST; 
hO2_ref_HGS=HGSsingle('H2','h',350,77.5);
hO2_inlet_HGS = hO2_ref_HGS - deltaH_O2; % kJ/mol

Hin = 2.8*hH2_inlet_HGS + 1.01*hO2_inlet_HGS;

species = {'H2', 'O2', 'H2O', 'OH', 'O', 'H'};
n0 = [2.8,1.01,0,0,0,0];
P0=62;

[Tp,~,np,~] = HGStp(species,n0,'H',Hin,P0);

%% The expansion
p = 60:-1:1;
[T_p_S,M_p_S,n_p_S,S_pre_S,S_post_S,H_pre_S,H_post_S,G_pre_S,G_post_S] = EXPANSION(species,Tp,P0,np,p,'Shifting');
[T_p_F,M_p_F,n_p_F,S_pre_F,S_post_F,H_pre_F,H_post_F,G_pre_F,G_post_F] = EXPANSION(species,Tp,P0,np,p,'Frozen');

%% Plots

figure(1)
plot(p,M_p_S,'b-')
hold on
plot(p,M_p_F,'r-')
grid
legend('Shifting','Frozen')
title('P vs Mach')

figure(2)
plot(p,T_p_S,'b-')
hold on
plot(p,T_p_F,'r-')
grid
legend('Shifting','Frozen')
title('P vs T')

n_sorted_F = cell(length(n_p_F{1}),1);
n_sorted_S = cell(length(n_p_F{1}),1);

for ii=1:length(n_p_F)
    for jj=1:length(n_p_F{ii})
       n_sorted_F{jj}(ii) =  n_p_F{ii}(jj);
       n_sorted_S{jj}(ii) =  n_p_S{ii}(jj);
    end
end

figure(3)
hold on
for ii=1:length(n_sorted_F)
    plot(p,n_sorted_F{ii},'-')
end
grid
legend(species)
title('n in Frozen')


figure(4)
hold on
for ii=1:length(n_sorted_S)
    plot(p,n_sorted_S{ii},'-')
end
grid
legend(species)
% ylim([5, 10])
title('n in Shifting')

figure(5)
plot(p,S_pre_S,'b-')
hold on
plot(p,S_pre_F,'r-')
plot(p,S_post_S,'g-')
plot(p,S_post_F,'c-')
grid
legend('S pre Shifting','S pre Frozen','S post Shifting','S post Frozen')
title('P vs S')

figure(6)
plot(p,S_pre_S.*T_p_S,'b-')
hold on
plot(p,S_pre_F.*T_p_F,'r-')
plot(p,S_post_S.*T_p_S,'g-')
plot(p,S_post_F.*T_p_F,'c-')
grid
legend('S pre Shifting','S pre Frozen','S post Shifting','S post Frozen')
title('P vs T*S')

figure(7)
plot(p,G_pre_S,'b-')
hold on
plot(p,G_pre_F,'r-')
plot(p,G_post_S,'g-')
plot(p,G_post_F,'c-')
grid
legend('G pre Shifting','G pre Frozen','G post Shifting','G post Frozen')
title('P vs G')

figure(8)
plot(p,H_pre_S,'b-')
hold on
plot(p,H_pre_F,'r-')
plot(p,H_post_S,'g-')
plot(p,H_post_F,'c-')
grid
legend('H pre Shifting','H pre Frozen','H post Shifting','H post Frozen')
title('P vs H')



%% The expansion function
function [T_p,M_p,n_p,S_pre,S_post,H_pre,H_post,G_pre,G_post] = EXPANSION(species,Tp,P0,np,p,Fro_Shift)
% To calculate the velocity we require the enthalpy per kilo so we convert
% and summ all mass of the inlet species.
% Inlet enthalpy per kilo and entropy
[S1,H1,MM] = HGSprop(species,np,Tp,P0,'S','H','Mm');
mt = sum(np)*MM/1000;
h1 = H1/mt;


if strcmp(Fro_Shift,'Shifting')
    exp_model = 1;
elseif strcmp(Fro_Shift,'Frozen')
    exp_model = 0;
else
    error('Ups,...The expansion model that you selected is not Frozen or Shifting. Check it pls')
end

% The points of P and starting points
n = np;
Tstar = Tp;
Pstar = P0;


for ii = 1:length(p)
    Pstar = p(ii);
    fprintf("P = %f\n",Pstar)
    
    if ii == 1
        opt = struct('xmin',300,'xmax',6000,'maxiter',200,'epsx',0.0001,'epsy', 0.0001,'fchange',1,'info',0,'dTp',50);
    else
        opt = struct('xmin',T_p(ii-1)-500,'xmax',T_p(ii-1)+100,'maxiter',200,'epsx',0.0001,'epsy', 0.0001,'fchange',1,'info',0,'dTp',50);
    end
    
    [Tstar,n] = HGSsecant(@hastobeS,n,opt);
    T_p(ii) = Tstar;
    [a,H2] = HGSprop(species,n,Tstar,Pstar,'a','H');   
    v2=sqrt(2*1000*(h1-H2/mt));
    n_p{ii} = n;
    M_p(ii) = v2/a;
    zeroM = 1 - v2/a;
%   disp(zeroM)
end

function [zeroS,n] = hastobeS(T,n)
%   fprintf("T = %f\n",T)
    [S_pre(ii),G_pre(ii),H_pre(ii)] = HGSprop(species,n,T,Pstar,'S','G','H'); 
    if exp_model
       [species,n,~] = HGSeq(species,n,T,Pstar);
    end

    [S_post(ii),G_post(ii),H_post(ii)] = HGSprop(species,n,T,Pstar,'S','G','H'); 
%   fprintf("S = %f\n\n",(S - S1))
    zeroS = (S_post(ii) - S1);   
end   
end

