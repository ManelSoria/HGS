%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC    
%***********************************************************************************************************
%
% Example 04: Adiabatic flame temperature with dissociation
%             this code is an example of how hgsTp works

function Ex04b_adiabatic_flame_T

clear; clc;

species={'H2','O2','H2O','H','O','OH'};
Tr=350;             % K 
P=10;               % bar
nr=[2;1;0;0;0;0];   % mol

Hin = HGSprop(species,nr,Tr,P,'H')


% Plot products enthalpy vs. T
T=linspace(300,5000,10);
Hout = zeros(length(T),1);
for i=1:length(T) 
    [~,comp,~] = HGSeq(species,nr,T(i),P);
    Hout(i) =HGSprop(species,comp,T(i),P,'H');    
end
plot(T,Hout,'-or',T,Hin*ones(length(T),1),'-ob');
legend('Hout','Hin'); 
xlabel('Temperature (K)'); ylabel('Enthalpy (kJ/molK)');

% Function to be solved to find T so that deltaH=0
    function DH=deltaH(Tp)
        [~,neq,~] = HGSeq(species,nr,Tp,P); % find equilibrium composition
        Ho = HGSprop(species,neq,Tp,P,'H');
        DH=Ho-Hin;
    end

% examples
fprintf('deltaH @ 400 K = %.2f kJ/molK \n',deltaH(400));
fprintf('deltaH @ 4000 K = %.2f kJ/molK \n',deltaH(4000));

% solving the equation
Tflame=fzero(@deltaH,3000,optimset('Display','iter'));

fprintf('Adiabatic flame temperature Tp = %.2f K \n',Tflame);

end
