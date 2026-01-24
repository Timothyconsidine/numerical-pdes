%Initializing bvp
domain = [-0.5 0.5];
kappa = @(x,u) 1 + 0*x;
f = @(x,u) -cos(u)+24.*x+cos(4.*x.^3-x);
d2f = @(x,u) sin(u);
bvp7 = TwoPointBVP(domain,kappa,true);
bvp7 = bvp7.set_f(f,d2f,false);
bvp7 = bvp7.set_left_boundary(2,-1,false);
bvp7 = bvp7.set_right_boundary(-2,-1,false);

%Initializing partition
s = 11;
N = 2^(s+3);
x = linspace(-0.5,0.5,N)';
partition = Partitions(x);

%Initializing FVEM System
fvem7 = FVEMSystem(partition,bvp7);
spatial7 = FVEMSystemMap(fvem7);

%Fixed-point iteration
maxnumiter = 100;
TOL = 1e-9;
initguess = zeros(N,1);
[alpha7, itercount, diff]  = spatial7.Iterate(initguess,maxnumiter,TOL);

%Plotting true vs approximate solution
y= linspace(-0.5, 0.5, 2^15+1);
truesol7 = 4.*y.^3-y;
clear figure(1);
hfig = figure(1);
fname = 'solplot7.pdf';
L = plot(y, truesol7, 'r', x, alpha7, 'b');
title("True vs. Approximate Solution for BVP 3")
xlabel('$x$');
ylabel('Solution');
legend('$u$','$u_h$, $N=2^{14}$');
 L(1).LineWidth = 3;
 L(2).LineWidth = 1;
plotsetting(hfig)
exportgraphics(figure(1),fname,'BackgroundColor','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error Calculation
u7 = @(x) 4.*x.^3-x;
Ns = @(s) 2^(s+3);
SE7 = ErrorCalculation3(bvp7,u7,Ns);
SN7 = [2^4 2^5 2^6 2^7 2^8 2^9 2^10 2^11 2^12 2^13 2^14];

%Error Plot
clear figure(2)
hfig2 = figure(2);
fname2 = 'logplot7.pdf';
loglog(SN7,SE7,'Color','red','Marker','*')
title("$\{N(s),E(s)\}: s = 1,...,11$ for BVP 3 (Log-Log Scale)")
xlabel('Partition Intervals: $N=2^{s+11}$')
ylabel("Error : $Norm(u-u_h)$")
plotsetting(hfig2)
exportgraphics(figure(2),fname2,'BackgroundColor','none');

%Regression
LM = fitlm(log(SN7),log(SE7), 'linear');
LM.Coefficients.Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%