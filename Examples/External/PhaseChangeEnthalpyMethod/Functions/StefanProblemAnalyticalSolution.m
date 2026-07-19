function [T,X] = StefanProblemAnalyticalSolution( T_L,T_S,T_m,L,rho,k_L,k_S,heatCapacity_L,heatCapacity_S,pos_x,time)
%StefanProblemAnalyticalSolution returns the analytical solution (i.e. 
% temperature and phase change Interface position) of the two-phase Stefan
% problem

% calculate the thermal diffusivity (liquid phase)
thermalDiffusivity_L=k_L/(rho*heatCapacity_L);

% calculate the thermal diffusivity (solid phase)
thermalDiffusivity_S=k_S/(rho*heatCapacity_S);

nu=sqrt(thermalDiffusivity_L/thermalDiffusivity_S);

% calculate the Stefan Number liquid
St_L=heatCapacity_L*(T_L-T_m)/L;

% calculate the Stefan Number solid
St_S=heatCapacity_S*(T_m-T_S)/L;

% calculate lambda by an approximation
trancendalEq = @(x) St_L/(exp(x^2)*erf(x))-St_S/(nu*exp((nu^2)*x^2)*erfc(nu*x))-x*sqrt(pi);
lambda=fzero(trancendalEq,0.3799,optimset('Display','notify'));

X_eq=@(t) 2*lambda*sqrt(thermalDiffusivity_L*t);
X=X_eq(time);

T=zeros(length(pos_x),1);
for i=1:length(pos_x)
    T(i)=calcT(T_L,T_S,T_m,pos_x(i),time,X,thermalDiffusivity_L,thermalDiffusivity_S,lambda);
end

end

function T=calcT(T_L,T_S,T_m,x,t,X,thermalDiffusivity_L,thermalDiffusivity_S,lambda)
    
    if x==X && x==0
        T=T_S;
    elseif x<X
        T=T_L-(T_L-T_m)*(erf(x/(2*sqrt(thermalDiffusivity_L*t))))/erf(lambda);
    else
        T=T_S+(T_m-T_S)*(erfc(x/(2*sqrt(thermalDiffusivity_S*t))))/erfc(lambda*sqrt(thermalDiffusivity_L/thermalDiffusivity_S));
    end

end