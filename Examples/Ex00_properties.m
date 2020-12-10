%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%* Commented by Paulino Gil -> Only this comment
%
%* ESEIAAT UPC          
%***********************************************************************************************************
%
% Example 00: Properties of single elements and mixtures

clear

Ru=8.3144621/1000; % kJ/molK


% Compute the propeties of a mixture of 1 mol of N2, 2 mol of CH4 and 3 mol
% of C3H8 at 3 bar and 400K 
[Mm,Cp,Cv,H,S,G,Rg,gamma,a]=HGSprop({'N2','CH4','C3H8'},[1,2,3],400,3); 
fprintf('propeties of a mixture of 1 mol of N2 , 2 mol of CH4 and 3 mol of C3H8 at 3 bar and 400K \n');
fprintf('Mm = %f, Cp = %f, Cv = %f \n',Mm,Cp,Cv);



% Using N2, verify that deltaH is aprox. equal to Cp*deltaT for small
% deltaT
[CpN2,hN2_1]=HGSprop({'N2'},1,300,10,'Cp','H');
[hN2_2]=HGSprop({'N2'},1,301,10,'H');

fprintf('Single species  N2 \n')
fprintf('DeltaH = <%.4f> , Cp = <%.4f> \n',(hN2_2-hN2_1),CpN2)




[Cpmix,hmix_1]=HGSprop({'N2','O2'},[1 1],300,10,'Cp','H');
[hmix_2]=HGSprop({'N2','O2'},[1 1],301,10,'H');

fprintf('Mixture of N2 and O2 \n')
fprintf('DeltaH = <%f> , Cp = <%f> \n',(hmix_2-hmix_1),Cpmix)




% Using N2, verify that S2-S1 aprox. = cp.ln(T2/T1)-Ru*ln(P2/P1)

[S1]=HGSprop({'N2'},1,400,20,'S');
[CpN2,RN2,S2]=HGSprop({'N2'},1,500,10,'Cp','Rg','S');

fprintf('Single species  N2 \n')
fprintf('S2-S1 = <%f> , cp*ln(T2/T1)-Ru*ln(P2/P1) = <%f> \n',S2-S1,CpN2*log(500/400)-Ru*log(10/20))




% Using air (80% N2, 20% O2) 
% verify @300K deltaH aprox. equal to Cp*deltaT
[CpA,hA_1]=HGSprop({'N2' 'O2'},[0.8 0.2],300,10,'Cp','H');
[hA_2]=HGSprop({'N2' 'O2'},[0.8 0.2],310,10,'H');

fprintf('Mixture of N2 and O2 \n')
fprintf('dH/10 = <%f> , Cp = <%f> \n',(hA_2-hA_1)/10,CpA)



% @300K and 1 bar verify:
T=300; % K
P=1; % bar
m=10; % kg (random)

[MM,Cp,Cv,H,S,G,Rg,gamma,a]=HGSprop({'N2' 'O2'},[0.8 0.2],T,P);

% sound speed  (J/kgK)^(1/2)=m/s
fprintf('Mixture of N2 and O2 \n')
fprintf('a = <%f> , sqrt(gamma*Rg*1000*T) = <%f> \n',a,sqrt(gamma*Rg*1000*T))


% Rg=Ru/MM  kJ/kg K
% The code is printing here the gas coefficient R
fprintf('1000*Ru/MM = <%f> , Rg = <%f> \n',1000*Ru/MM ,Rg)


% Cp/Cv=gamma
fprintf('Cp/Cv-gamma = <%f> \n',Cp/Cv-gamma)


% g=h-T*s
v=1000*Rg*T/(P*1e5); % m^3/kg

V=m*v;
n=m*1000/MM; % mol

h=H/m; % kJ/kg
s=S/m; % kJ/kgK
g=G/m; % kJ/kg
u=h-P*1e5*v/1000;

fprintf('h-T*s = <%f>,g = <%f> \n',h-T*s,g)

% verify molar cp of the mixture
CpO2=HGSsingle({'O2'},'cp',T,1); % kJ/molK
MMO2=HGSsingle({'O2'},'Mm',T,1); 
CpN2=HGSsingle({'N2'},'cp',T,1); % kJ/molK
MMN2=HGSsingle({'N2'},'Mm',T,1); 
hN2=HGSsingle({'N2'},'h',T,1); % kJ/mol

fprintf('Verifying independently with HGSsingle (0.8*CpN2+0.2*CpO2-Cp) = <%f> \n',0.8*CpN2+0.2*CpO2-Cp)


CpN2kg=1000*CpN2/MMN2; % kJ/molK -> kJ/kgK
hN2kg=1000*hN2/MMN2; % kJ/mol -> kJ/kgK 

fprintf('N2: Cp=%f kJ/kgK hN2 =%f kJ/kg \n',CpN2kg, hN2kg );


% Note that the program returns in units of kJ/molK. Here are converted to
% kJ/kgK
Cpkg=Cp/(MM/1000);

fprintf('cp=%f kJ/kgK \n',Cpkg);

% hgssingle computes H, G and S properties for a single element
fprintf('hgssingle computes H, G and S properties for a single element \n')
HGSsingle('O2','h',500,10)
 
% hgsprop is also able to compute for a mixture of gases. The sintax is
% shown below
fprintf('hgsprop is also able to compute for a mixture of gases \n')
[Mm,Cp,Cv,H,S,G,Rg,gamma,a]=HGSprop({'O2' 'N2'},[0.2 0.8],400,3)

fprintf('Mm=%f \n', Mm);