% data1: sea-water Pc
% Written by Ali A. Eftekhari
% Last checked: June 2021
data{1}=[0.073 249
0.085 45
0.09 15
0.131 8
0.197 3.33
0.23 2.22
0.25 1.79
0.328 0.6
0.38 0
0.5 -0.1
0.6 -0.3
0.7 -0.7
0.8 -2.5
0.831 -5.8
0.84 -10
0.856 -20
0.863 -32.2
0.867 -43
0.869 -49.1
0.87 -54.3
0.882 -92.8
0.893 -126.9
0.895 -142.5
0.899 -167.7
0.901 -250];
% zero sulphate Pc, table 3
data{2}=[0.056 200
0.058 36.9
0.073 9.3
0.112 5.06
0.137 3.56
0.176 1.46
0.232 0.5
0.28 0
0.503 -0.8
0.606 -1.13
0.698 -2
0.768 -3
0.793 -4
0.81 -5
0.827 -13
0.84 -31.4
0.857 -79.3];
% JCR library for cretaceous chalk
data{3}=[0.079 0.1
0.1 0.02
0.3 0.015
0.54 0.0
0.68 -0.4
0.747 -0.7];
% fit model to data
pc1=data{3}(:,2)*1e5; % Pa
sw1=data{3}(:,1);
sw_plot=linspace(min(sw1),max(sw1),5000);
% pc_plot=interp1(sw, pc, sw_plot, 'linear', 'extrap');
% plot(sw, pc, 'o', sw_plot, pc_plot)
% plot(diff(pc_plot)./diff(sw_plot))
% log(pc)=log(pce)-(1/labda)*log((sw-sw0)/(1-sw0-sor))
% ind0=find(pc==0, 1);
% p1=polyfit(sw(1:ind0-1), log(pc(1:ind0-1)), 1)
% plot(sw(1:ind0-1), log(pc(1:ind0-1)), 'o')

pc=@(sw, cwi, coi, awi, aoi, swc, sor)(cwi./((sw-swc)/(1-swc)).^awi+coi./((1-sw-sor)/(1-sor)).^aoi);
swc0=0.0789;
sor0=0.2529;
labda=2.4;
f=@(x, sw)pc(sw, x(1), x(2), 1/labda, 1/labda, swc0, sor0);
fw=@(x, sw)pc(sw, x(1), 0.0 , x(3), x(4), swc0, sor0);
fo=@(x, sw)pc(sw, 0.0, x(2), x(3), x(4), swc0, sor0);
x=lsqcurvefit(f, [1e3, -1e3],sw1, pc1)
% fminsearch (base MATLAB) is used instead of patternsearch so this example
% does not require the Global Optimization Toolbox.
x=fminsearch(@(x)(sum(abs(f(x, sw1)-pc1))), x)
plot(sw1, pc1, 'o', sw_plot, f(x, sw_plot))%, ...
% sw_plot, fw(x, sw_plot), '--', ...
% sw_plot, fo(x, sw_plot), '-.'); 
grid;