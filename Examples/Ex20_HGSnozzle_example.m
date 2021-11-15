% HGSnozzle example

clear
close all
sp = {'H2', 'O2', 'H2O', 'OH', 'O', 'H'};
Pchamber=62; % bar
Pe=0.4; % bar
Pa=1; % bar
Tinlet=300; % K

n0=[2 1 0 0 0 0 ]*1e1; % inlet to combustion chamber
mfr0=sum(n0)*HGSprop(sp,n0,Tinlet,Pchamber,'Mm')*1e-3; % kg/s = (mol/s)* (kg/mol)

[Tcomb,sp,n] = HGStp(sp,n0,'T',Tinlet,Pchamber)

% nozzle inlet mass flow rate

mfr1=sum(n)*HGSprop(sp,n,Tcomb,Pchamber,'Mm')*1e-3; % kg/s = (mol/s)* (kg/mol)


NP=20;
PV=linspace(Pchamber,Pe,NP);
NV=zeros(NP,numel(sp)); % NV(point,element)
for i=1:numel(PV)
    fprintf('Pressure %d / %d \n',i,NP);
    [TV(i),sp,NV(i,:),M(i),flag]=HGSisentropic(sp,n,Tcomb,Pchamber,'Shifting','P',PV(i));
    % check flag
    [MM(i),Rg(i),a(i),H(i)]=HGSprop(sp,NV(i,:),TV(i),PV(i),'Mm','Rg','a','H');
    v(i)=M(i)*a(i); % velocity m/s
    rho(i) = PV(i)/(Rg(i)*TV(i));
    NT(i)=sum(NV(i,:)); % total number of mol at point i
    mfr(i)=NT(i)*MM(i)*1e-3; % mass flow rate at point i (has to be constant)
    A(i)=mfr(i)/(rho(i)*v(i));
    F(i)=mfr(i)*v(i)+A(i)*(PV(i)-Pa)*1e5; % N UHHHH !! 
end

%%
close all

plot(Pchamber-PV,TV);
title('Temperature vs P');

figure

plot(Pchamber-PV,rho);
title('rho vs P');

figure
plot(Pchamber-PV,NV(:,1)./NT(:)); % H2 fraction
title('H2 fraction vs P');

figure
plot(Pchamber-PV,v); 
title('Velocity vs P');


figure
plot(Pchamber-PV,M); 
title('M vs P');

figure
plot(Pchamber-PV,A); 
title('A vs P');

figure
plot(Pchamber-PV,F); 
title('F vs P');