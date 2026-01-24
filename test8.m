%Initializing bvp
domain = [0 1];
kappa = @(x,u) 2 + 2.*u;
f = @(x,u) -3.*x.^2;
d2f = @(x,u) 0*x;
bvp8 = TwoPointBVP(domain,kappa,false);
bvp8 = bvp8.set_f(f,d2f,true);
bvp8 = bvp8.set_left_boundary(0,0,true);
bvp8 = bvp8.set_right_boundary(0,0,true);

%Initializing partition
s = 11;
N = 2^(s+3);
x = linspace(0,1,N)';
partition = Partitions(x);

%Initializing FVEM System
fvem8 = FVEMSystem(partition,bvp8);
spatial8 = FVEMSystemMap(fvem8);

%Fixed-point iteration
initguess = zeros(N,1);
maxnumiter = 100;
TOL = 1e-9;
[alpha8, itercount, diff]  = spatial8.Iterate(initguess,maxnumiter,TOL);

%Plotting true vs approximate solution
y= linspace(0, 1, 2^15+1);
truesol8 = @(y) -1 +1/2*sqrt(1./y.^2).*y.*sqrt(-y.^4 + y + 4);
clear figure(1);
hfig = figure(1);
fname = 'solplot8.pdf';
L = plot(y, truesol8(y), 'red', x, alpha8, 'blue');
title("True vs. Approximate Solution for BVP 4")
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^{14}$');
 L(1).LineWidth = 3;
 L(2).LineWidth = 1;
plotsetting(hfig)
exportgraphics(figure(1),fname,'BackgroundColor','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error Calculation
u8 = @(x) -1 +1/2*sqrt(1./x.^2).*x.*sqrt(-x.^4 + x + 4);
Ns = @(s) 2^(s+3);
SE8 = ErrorCalculation3(bvp8,u8,Ns)
SN8 = [2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 2^14];

%Error Plot
clear figure(2)
hfig2 = figure(2);
fname2 = 'logplot8.pdf';
loglog(SN8,SE8,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 4 (Log-Log Scale)")
xlabel('Partition Intervals: $N=2^{s+11}$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig2)
exportgraphics(figure(2),fname2,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN8),log(SE8), 'linear');
LM.Coefficients.Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%