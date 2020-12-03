%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC         
%***********************************************************************************************************
%
% Example 09: Comparison between hgs / rpa using RP-1
%             RP-1 // O2 isentropic expansion

clear; clc;

species={'CO','CO2','COOH','H','H2','H2O','H2O2','HCO-','HO2','O','O2','OH'};

% RPA has HCO, but not HCO-

n1=[0.2292308;...
    0.1658024;...
    0.0000007;...
    0.0447093;...
    0.0488284;...
    0.2672289;...
    0.0000017;...
    0.0000007;...
    0.0000481;...
    0.0477951;...
    0.1029129;...
    0.0934408
];

p1=1;
T1=3082.3911;

[~,neq,~] = HGSeq(species,n1,T1,p1)

norm(neq-n1)/norm(neq)

p2=0.1

[T2,n2,~,~,~,~] = HGSisentropic(species,n1,T1,p1,p2)

