%Initializing bvp
domain = [-1 1];
kappa = @(x,u) 1 + 0*x;
f = @(x,u) sin(2.*x).*(cos(x.^2-1)-cos(u))+2.*(u+2-x.^2);
d2f = @(x,u) sin(2.*x).*sin(u)+2;
bvp6 = TwoPointBVP(domain,kappa,true);
bvp6 = bvp6.set_f(f,d2f,false);
bvp6 = bvp6.set_left_boundary(0,0,true);
bvp6 = bvp6.set_right_boundary(0,0,true);

%Initializing partition
s = 11;
N = 2^(s+3);
x = linspace(-1,1,N)';
partition = Partitions(x);

%Initializing FVEM System
fvem6 = FVEMSystem(partition,bvp6);
spatial6 = FVEMSystemMap(fvem6);

%Fixed-point iteration
initguess = zeros(N,1);
maxnumiter = 100;
TOL = 1e-9;
[alpha6, itercount, diff]  = spatial6.Iterate(initguess,maxnumiter,TOL);

%Plotting true vs approximate solution
y= linspace(-1, 1, 2^15+1);
truesol6 = y.^2-1;
clear figure(1);
hfig = figure(1);
fname = 'solplot6.pdf';
L = plot(y, truesol6, 'r', x, alpha6, 'b');
title("True vs. Approximate Solution for BVP 2")
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^{14}$');
 L(1).LineWidth = 3;
 L(2).LineWidth = 1;
plotsetting(hfig)
exportgraphics(figure(1),fname,'BackgroundColor','none');

%Error Calculation
u6 = @(x) x.^2-1;
Ns = @(s) 2^(s+3);
SE6 = ErrorCalculation3(bvp6,u6,Ns);
SN6 = [2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 2^14];

%Error Plot
clear figure(2)
hfig2 = figure(2);
fname2 = 'logplot6.pdf';
loglog(SN6,SE6,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 2 (Log-Log Scale)")
xlabel('Partition Intervals: $N=2^{s+3}$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig2)
exportgraphics(figure(2),fname2,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN6),log(SE6), 'linear');
LM.Coefficients.Estimate