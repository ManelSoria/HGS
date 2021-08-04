function [Tp,n,flag] = HGSsecant(f,n0,options)
%**************************************************************************
%
% [Tp,n,flag] = HGSsecant(f,n0,options)
%
%**************************************************************************
%
% HGSsecant solves a function using a combination of the secant method and
% bisection method
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% f --> Function 
% n0 --> [mol] Initial mixture
% options --> Structure with the options for the secant method. 
%                 .xmin [K] Temperature minimum for the solver;
%                 .xmax [K] Temperature maximum for the solver;
%                 .maxiter Max iterations for the solver;
%                 .epsx Diferential T where the solver reachs the solution;
%                 .epsy Diferential S where the solver reachs the solution;
%                 .fchange T difference where secant method is
%                          changed by bisection method;
%                 .type Select between: 'Frozen' for frozen flow
%                                       'Shifting' for shifting flow
%                 .info Detailed info == 1; No info == 0.
%                 .dTp Improve the velocity with the approximation of
%                 parabola. +- dTp
%           struct('xmin',300,'xmax',6000,'maxiter',200,'epsx',0.1,'epsy',
%                   1,'fchange',5,'type','Shifting','info',0,'dTp',100)
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Final temperature
% n --> [mol] Final mixture
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

%% Options

def.xmin = 300;
def.xmax = 4000;
def.maxiter = 200;
def.epsx = 0.1;
def.epsy = 1;
def.fchange = 5;
def.info = 0;
def.dTp = 100;

if exist('options','var') && ~isempty(options)
    fields = fieldnames(options);
    for ii = 1:length(fields)
        def.(fields{ii}) = options.(fields{ii});
    end
end


x1 = def.xmin ;
x2 = def.xmax ;
maxiter = def.maxiter;
epsx = def.epsx;
epsy = def.epsy;
fchange = def.fchange;
info = def.info;
dTp = def.dTp;

%% Evaluation of Max and Min T

[y1,n1] = f(x1,n0); 
[y2,n2] = f(x2,n0); 

Tp=[];
n=[];

if y1*y2 > 0 % No sign change, sorry !
    flag=-2; % Initial sign change not found
    return
end

if x2-x1 > 1500  % Try to fit to a parabola and solve the eq.
   x3 = (x2+x1)/2;
   [y3,~] = f(x3,n0); 
   a = (y1 -(y2-y3)/(x2-x3)*x1 - y3 + x3*(y2-y3)/(x2-x3)) / (x1^2 + (x1-x3)*(x3^2-x2^2)/(x2-x3) - x3^2);
   b = (y2-y3+a*x3^2-a*x2^2) / (x2-x3);
   c = y3 - a*x3^2 - b*x3;
   
   % Solution of the parabola
   xsol(1) = (-b+sqrt(b^2-4*a*c))/(2*a);
   xsol(2) = (-b-sqrt(b^2-4*a*c))/(2*a);
   
   if xsol(1) >= x1
       x1p = xsol(1)-dTp;
       x2p = xsol(1)+dTp;
   else
       x1p = xsol(2)-dTp;
       x2p = xsol(2)+dTp;
   end
   
   % Don't allow the current range to extend the range specified
   if x1p < options.xmin
      x1p = options.xmin;
   end
   if x2p > options.xmax
      x2p = options.xmax;
   end
   
   [y1p,n1] = f(x1p,n0); 
   [y2p,n2] = f(x2p,n0); 
   
   if y1p*y2p < 0 % Parabola method is okey
        x1 = x1p;
        y1 = y1p;
        x2 = x2p;
        y2 = y2p;
   end
end

flag=-1; % We assume we are not solving it

for ii=1:maxiter
        
    if x2-x1 < fchange % If limits are close, switch to bisection algorithm
        xc = (x1+x2)/2; 
    else
        xc = x1-y1*(x2-x1)/(y2-y1); % Secant method
    end
   
    if xc-x1 < x2-xc % Next iteration is closer to x1
        n=n1;
    else
        n=n2;
    end
    
    [yc,n]=f(xc,n); % Compute next value
    
    if abs(yc)<epsy || (abs(xc-x1)<epsx && abs(x2-xc)<epsx )% Stop if it is solved
        flag = 1;
        Tp=xc;
        break;
    end
    
    if info
       fprintf('ii=%d x1=%e y1=%e xc=%e yc=%e x2=%e y2=%e \n',ii,x1,y1,xc,yc,x2,y2);
    end
    
    if yc*y1>0 % Change limits
        y1=yc;
        x1=xc;
        n1=n;
    else
        y2=yc;
        x2=xc;
        n2=n;
    end
    
end

end

