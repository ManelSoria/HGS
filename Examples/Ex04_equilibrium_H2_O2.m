%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC    
%***********************************************************************************************************
%
% Example 04: H2O equilibrium dissociation for different values of
%             temperature (T)
%
% H20 <-> H2 + O2 + H + O + OH

clear; clc;

p=1;                        % bar
T=linspace(300,5000,50);    % K

% loop to compute composition using hgseq
lenT = length(T);
nH2 = zeros(lenT,1);
nO2 = zeros(lenT,1);
nH2O = zeros(lenT,1);
nH = zeros(lenT,1);
nO = zeros(lenT,1);
nOH = zeros(lenT,1);
for i=1:lenT
    fprintf('Solving equilibrium composition for T=%f K\n',T(i));
    [~,comp,~]=HGSeq({'H2','O2','H2O','H','O','OH'},[2;1;0;0;0;0],T(i),1);
    nH2(i)=comp(1);
    nO2(i)=comp(2);
    nH2O(i)=comp(3);
    nH(i)=comp(4);
    nO(i)=comp(5);
    nOH(i)=comp(6);
end

% plot the results and apply some basic formatting
plot(T,nH2O,'r',T,nH,'b',T,nO,'g',T,nOH,'k',T,nH2,'c',T,nO2,'m');
legend('H20','H','O','OH','H2','O2','Location','NorthWest');
xlabel('Temperature (K)'); ylabel('concentration (mol)');
grid;
