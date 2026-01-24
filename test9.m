%Initializing bvp
domain = [1 6];
kappa = @(x,u) 1./(x.*exp(u));
f = @(x,u) -1+0*x;
d2f = @(x,u) 0*x;
bvp9 = TwoPointBVP(domain,kappa,false);
bvp9 = bvp9.set_f(f,d2f,true);
bvp9 = bvp9.set_left_boundary(0,0,true);
bvp9 = bvp9.set_right_boundary(6,0,false);

%Initializing partition
s = 11;
N = 2^(s+3);
x = linspace(1,6,N)';
partition = Partitions(x);

%Initializing FVEM System
fvem9 = FVEMSystem(partition,bvp9);
spatial9 = FVEMSystemMap(fvem9);

%Fixed-point iteration
initguess = ones(length(x),1).*-2.9;
maxnumiter = 100;
TOL = 1e-9;
[alpha9, itercount, diff]  = spatial9.Iterate(initguess,maxnumiter,TOL);

%Plotting true vs approximate solution
y= linspace(1, 6, 2^15+1);
truesol9 = log(3)-log(y.^3+2);
clear figure(1);
hfig = figure(1);
fname = 'solplot9.pdf';
L = plot(y, truesol9, 'r', x, alpha9, 'b');
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^s+3$');
 L(1).LineWidth = 3;
 L(2).LineWidth = 1;
plotsetting(hfig)
exportgraphics(figure(1),fname,'BackgroundColor','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error Calculation
u9 = @(x) log(3)-log(x.^3+2);
Ns = @(s) 2^(s+3);
SE9 = ErrorCalculation3(bvp9,u9,Ns);
SN9 = [2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 2^14];

%Error Plot
clear figure(2)
hfig2 = figure(2);
fname2 = 'logplot9.pdf';
loglog(SN8,SE8,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 4 (Log-Log Scale)")
xlabel('Partition Intervals: $N=2^{s+11}$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig2)
exportgraphics(figure(2),fname2,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN9),log(SE9), 'linear');
LM.Coefficients.Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%