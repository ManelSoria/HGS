function [A,b,Aeq,beq,lb,ub] = HGSparamFmincon(id,n0)
%**************************************************************************
%
% [A,b,Aeq,beq,lb,ub] = HGSparamFmincon(id,n0)
%
%**************************************************************************
%
% HGSparamFmincon returns the parameters for fmincon solver
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% id --> Id of species
% n0 --> [mol] Species mols
%
% Outputs:
%--------------------------------------------------------------------------
% A --> fmincon A inequality matrix A*x < b
% b --> fmincon b inequality vector A*x < b
% A --> fmincon Aeq equality matrix Aeq*x = beq
% b --> fmincon beq equality vector Aeq*x = beq
% lb --> fmincon lb lower boundary limit
% ub --> fmincon ub upper boundary limit
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Mir�
% *ESEIAAT UPC    

spec = length(id);
global HGSdata;HGSload

%% Inequality: Ax=<b.

A = [];
b = [];


%% Equality: Aeq*x=beq.
% Mass conservation of elemts must be complished
Elem = {};
for ii=1:spec
    for jj = 1:4
        Elem = [Elem; HGSdata.spec{id(ii),jj}];
    end
end

Elem = unique(Elem);
Lelem = length(Elem);
Aeq = zeros(Lelem,spec);
beq = zeros(Lelem,1);


for ii=1:Lelem
    for jj=1:spec
        for kk=1:4
            if strcmp(HGSdata.spec{id(jj),kk},Elem{ii,1})
                beq(ii,1) = beq(ii,1) + n0(jj)*HGSdata.nspec{id(jj),kk};
                Aeq(ii,jj) = HGSdata.nspec{id(jj),kk};
            end
        end
    end
end




%% Lower limit: lb && Upper limit: ub

lb = zeros(spec,1);
ub = [];


end