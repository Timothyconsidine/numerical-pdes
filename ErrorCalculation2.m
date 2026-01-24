function ES = ErrorCalculation2(bvp, truesol, N)
    
    E = zeros(11,1); %Initializing vector of error values

    for s = 1:11 %Calculating the error for each partition size

        K = N(s);

        %Constructing non-uniform partition
        a = linspace(-1,.07,2000*K+1);
        b = linspace(.07,.99,8*K+1);
        c = linspace(.99,1,40*K+1);
        x = [a,b(2:end),c(2:end)]';
        partition = Partitions(x);

        initguess = sparse((2000*K+1)+(8*K)+(40*K),1);
        fvem = FVEMSystem(partition,bvp);
        fvemmap = FVEMSystemMap(fvem);
        
        %Solving for approximate solution
        u_h = fvemmap.Iterate(initguess,100,1e-7); 
        x = x';
        u = truesol(x)';
        
        %Caclulating Error
        E(s) = sqrt((1/2)*sum(((u(2:length(u_h))-u_h(2:length(u_h))) ...
            -(u(1:length(u_h)-1)-u_h(1:length(u_h)-1))).^2));

    end
    
    ES = E; %Returning vector of error values
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%