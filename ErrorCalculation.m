function ES = ErrorCalculation(bvp, truesol, N)
    
    E = zeros(11,1); %Initializing vector of error values

    for s = 1:11 %Calculating the error for each partition size
        
        K = N(s);
        x = linspace(bvp.domain(1),bvp.domain(2),K)'; 
        initguess = zeros(K,1);
        
        %Constructing partition
        partition = Partitions(x); 
        fvem = FVEMSystem(partition,bvp);
        fvemmap = FVEMSystemMap(fvem);

        %Solving for approximate solution
        u_h = fvemmap.Iterate(initguess,100,1e-9); 
        x = x';
        u = truesol(x)';
        %Caclulating Error
        E(s) = sqrt((1/2)*sum(((u(2:length(u_h))-u_h(2:length(u_h))) ...
            -(u(1:length(u_h)-1)-u_h(1:length(u_h)-1))).^2));

    end
    
    ES = E; %Returning vector of error values
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%