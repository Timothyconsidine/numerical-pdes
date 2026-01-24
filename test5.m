%Initializing bvp
domain = [-1 1];
kappa = @(x,u) 10^(-6)+0*x;
f = @(x,u) x.*u;
d2f = @(x,u) x;
ga = 1;
gb = 1;
bvp5 = TwoPointBVP(domain,kappa,true);
bvp5 = bvp5.set_f(f,d2f,true);
bvp5 = bvp5.set_left_boundary(ga,0,true);
bvp5 = bvp5.set_right_boundary(gb,0,true);

%Initializing Partition
s = 11;
x = linspace(-1,1,2^(s+11))';
partition = Partitions(x);

%Initializing FVEM System
fvem5 = FVEMSystem(partition,bvp5);
spatial5 = FVEMSystemMap(fvem5);

%Fixed-point iteration
initguess = zeros(2^(s+11),1);
maxnumiter = 100;
TOL = 1e-7;
[alpha5, itercount, diff]  = spatial5.Iterate(initguess,maxnumiter,TOL);

%True Solution
A(1,:) = [airy(-100) airy(2,-100)];
A(2,:) = [airy(100) airy(2,100)];
B = [1 1]';
c = A\B;

%Plotting True vs Approximate Solution for uniform partition
y = linspace(-1, 1, 2^22+1);
truesol5 = c(1)*airy(100*y)+c(2)*airy(2,100*y);
clear figure(1);
hfig = figure(1);
fname = "solplot5a.pdf";
L = plot(y, truesol5, 'r', x, alpha5, 'b');
title("True vs. Approximate Solution for BVP 1")
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^{22}$');
 L(1).LineWidth = 3;
 L(2).LineWidth = 1;
plotsetting(hfig)
exportgraphics(figure(1),fname,'BackgroundColor','none');

%Error Calculation for uniform partition
u5 = @(x) c(1)*airy(100*x)+c(2)*airy(2,100*x);
Ns = @(s) 2^(s+11);
SE5 = ErrorCalculation(bvp5,u5,Ns);
SN5 = [2^12 2^13 2^14 2^15 2^16 2^17 2^18 2^19 2^20 2^21 2^22];

%Error plot for uniform partition
clear figure(2)
hfig2 = figure(2);
fname2 = 'logplot5a.pdf';
loglog(SN5,SE5,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 1 (Log-Log Scale)")
xlabel('Partition Intervals: $N=2^{s+11}$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig2)
exportgraphics(figure(2),fname2,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN5),log(SE5), 'linear');
LM.Coefficients.Estimate

%Creating non-uniform partition
a = linspace(-1,.07,2000*2^s+1);
b = linspace(.07,.99,8*2^s+1);
d = linspace(.99,1,40*2^s+1);
k = [a,b(2:end),d(2:end)]';
partitiona = Partitions(k);

%Initializing FVEM System using the non-uniform partition
fvem5a = FVEMSystem(partitiona,bvp5);
spatial5a = FVEMSystemMap(fvem5a);

%Fixed-point Iteration
initguess2 = zeros((2000*(2^s)+1)+(8*(2^s))+(40*(2^s)),1);
[alpha5a, itercount2, diff2]  = spatial5a.Iterate(initguess2,maxnumiter,TOL);
alpha5a = alpha5a';

%Plotting True vs Approximate Solution for non-uniform partition
p = linspace(-1,1,(2000*(2^s)+1)+(8*(2^s))+(40*(2^s)));
truesol5a = c(1)*airy(100*p)+c(2)*airy(2,100*p);
clear figure(3);
hfig3 = figure(3);
fname3 = "solplot5b.pdf";
T = plot(p, truesol5a, 'r',k, alpha5a, 'b');
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^{s+11}$');
 T(1).LineWidth = 3;
 T(2).LineWidth = 1;
plotsetting(hfig3)
exportgraphics(figure(3),fname3,'BackgroundColor','none');

%Error Calculation for non-uniform partition
Nsa = @(s) 2^s;
SE5a = ErrorCalculation2(bvp5,u5,Nsa);
SN5a = [2^1 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11];

%Error plot for non-uniform partition
clear figure(4)
hfig4 = figure(4);
fname4 = 'logplot5b.pdf';
loglog(SN5a,SE5a,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 1 (Log-Log Scale)")
xlabel('Partition Intervals: $N=\sum^{4}_{j=1}N_j(s)$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig4)
exportgraphics(figure(4),fname4,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN5a),log(SE5a), 'linear');
LM.Coefficients.Estimate
