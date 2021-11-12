%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    
%***********************************************************************************************************
%
% Adiabatic H2 / O2 reaction with dissociation
% (as an example of how HPStp looks for the temperature)

% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear
close 

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=66; % bar
nr=[2;1;0;0;0;0]; % mol
 
Hin = HGSprop(species,nr,Tr,P,'H');
 
T=linspace(300,5000,20);
for i=1:length(T) 
    [~,comp,~]= HGSeq(species,nr,T(i),P);
    Hout(i) = HGSprop(species,comp,T(i),P,'H');    
end
 
plot(T,Hout,'-r',T,Hin*ones(length(T),1),'-b','LineWidth',2)
hold on


[Tp,~,np,~] = HGStp(species,nr,'T',Tr,P);
flameHout = HGSprop(species,np,Tp,P,'H');
plot(Tp,flameHout,'ok','MarkerSize',14,'MarkerFaceColor','k');

set(gca,'FontSize',18)
xlabel('T');
ylabel('H');
legend({'Hout','Hin','HGStp reaction temperature'});
title('H2 O2 combustion with dissociation');