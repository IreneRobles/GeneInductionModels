# Function that generates the transition probability matrix
# Very prout to say that today, 26 Oct 2018, the matrix seem to work

function A_matrix_feedback_0mRNASTART(parameters, N)
    
    k_a = parameters[1]
    k_i = parameters[2]
    k_on = parameters[3]
    k_off = parameters[4]
    k_tc = parameters[5]
    k_dg = parameters[6]
    
    
    # 4 possible states 0, 1, or , 2 alleles activated, or inactive state
    total_states = 4 
    
    #### This is for simple telegraph model 
    ## A = zeros(Float64, total_states*N, total_states*N)
    #### In this case we need an extra state which we are going to put at the end
    
    A = zeros(Float64, total_states*N+1, total_states*N+1)

    ## Assemble A row by row
    
    #####################################################################################
    ### Special case is when the cell has 0 molecules of mRNA (it cannot degrade any) ###
    #####################################################################################
    ### For CI0
    A[end,end] = -k_a 
    
    ### For CA0
    # Cell goes CA0(0) --> CA1(0) 2k_on since cell is doploid
    A[1,1] = -2*k_on - k_i #
    # Cell goes CA0(1) --> CA0(n) 
    A[1,2] = k_dg #
    # Cell goes CA1(0) --> CA0(0)
    A[1,1+N] = k_off #
    # Cell goes CI0(0) --> CA(0)
    A[1,end] = k_a 
    
    ### For CA1
    # Cell goes CA0(n) --> CA1(n)
    A[1+N,1] = 2*k_on #
    # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n)
    A[1+N,1+N] = -(k_off + k_tc + k_on + k_i) #
    # Cell goes CA1(n + 1) --> CA1(n)
    A[1+N,2+N] = k_dg #
    # Cell goes CA2(n) --> CA1(n)
    A[1+N,1+2N] = 2*k_off #
    
    ### For CA2
    # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
    A[1+2N,1+2N] = -(2*k_tc + 2*k_off + k_i)
    # Cell goes CA2(n + 1) --> CA2(n)
    A[1+2N,2+2N] = k_dg #
    # Cell goes CA1(n) --> CA2(n)
    A[1+2N,1+N] = k_on #
    
    ### For CJ0
    # Cell goes CJ0(n+1) --> CJ0(n) 
    A[1+3N,2+3N] = k_dg
    # Cell goes CA0(n) --> CJ0(n)
    A[1+3N,1] = k_i 
    # Cell goes CA1(n) --> CJ0(n))
    A[1+3N,1+N] = k_i 
    # Cell goes CA2(n) --> CJ0(n))
    A[1+3N,1+2N] = k_i 
    
    #####################################################################################

    for ii = 2:N-1
        ### For CA0
        # Cell goes CA0(n) --> CA1(n) or CA0(n) --> CA0(n - 1) 
        # ii - 1 because the ii correspond to 0 molecules
        A[ii,ii] = -(2*k_on+k_dg*(ii-1)+ k_i)
        # Cell goes CA0(n+1) --> CA0(n)
        A[ii,ii+1] = k_dg*ii #
        # Cell goes CA1(n) --> CA0(n)
        A[ii,ii+N] = k_off #
        
        ### For CA1
        # Cell goes CA0(n) --> CA1(n)
        A[ii+N,ii] = 2*k_on #
        # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n) or CA1(n) --> CA1(n - 1)
        A[ii+N,ii+N] = -(k_tc+k_off+k_on+k_dg*(ii-1) + k_i)
        # Cell goes CA1(n+1) --> CA1(n)
        A[ii+N,ii+N+1] = k_dg*ii # 
        # Cell goes CA1(n-1) --> CA1(n)
        A[ii+N,ii+N-1] = k_tc #
        # Cell goes CA2(n) --> CA1(n)
        A[ii+N, ii+2N] = 2*k_off #
        
        
        ### For CA2
        # Cell goes CA1(n) --> CA2(n)
        A[ii+2N,ii+N] = k_on #
        # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1) or CA2(n) --> CJ0(n)
        A[ii+2N,ii+2N] = -(2*k_tc+2*k_off+k_dg*(ii-1) + k_i) #
        # Cell goes CA2(n+1) --> CA2(n)
        A[ii+2N,ii+2N+1] = k_dg*ii #
        # Cell goes CA2(n-1) --> CA2(n)
        A[ii+2N,ii+2N-1] = 2*k_tc #
        
        ### For CJ0
        # Cell goes CJ0(n+1) --> CJ0(n) 
        A[ii+3N,ii+3N + 1] = k_dg*ii
        # Cell goes CJ0(n) --> CJ0(n+1) 
        A[ii+3N,ii+3N] = -k_dg*(ii-1)
        # Cell goes CA0(n) --> CJ0(n)
        A[ii+3N,ii] = k_i 
        # Cell goes CA1(n) --> CJ0(n))
        A[ii+3N,ii+N] = k_i 
        # Cell goes CA2(n) --> CJ0(n))
        A[ii+3N,ii+2N] = k_i 
        
         
    end
    
    
    #####################################################################################
    ### Special case is when the cell has maximum molecules of mRNA (it cannot degrade any) ###
    #####################################################################################
   
    ### For CA0
    # Cell goes CA0(n) --> CA1(n) or CA0(n) --> CA0(n - 1) 
    A[N,N] = -(2*k_on+k_dg*(N-1) + k_i)
    # Cell goes CA1(n) --> CA0(n)
    A[N,2*N] = k_off
    
            
    ### For CA1
    # Cell goes CA0(n) --> CA1(n)
    A[2*N,N] = 2*k_on
    # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n) or CA1(n) --> CA1(n - 1)
    A[2*N,2*N] = -(k_off+k_on+k_dg*(N-1)+k_i)
    # Cell goes CA1(n-1) --> CA1(n)
    A[2*N,2*N-1] = k_tc
    # Cell goes CA2(n) --> CA1(n)
    A[2*N, 3*N] = 2*k_off
    
    ### For CA2
    # Cell goes CA1(n) --> CA2(n)
    A[3*N,2*N] = k_on
    # Cell goes  CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
    A[3*N,3*N] = -(2*k_off+k_dg*(N-1)+k_i)
    # Cell goes CA2(n-1) --> CA2(n)
    A[3*N,3*N-1] = 2*k_tc
    
    ### For CJ0
    # Cell goes CJ0(n) --> CJ0(n+1) 
    A[4N,4N] = -k_dg*(N-1)
    # Cell goes CA0(n) --> CJ0(n)
    A[4N,N] = k_i 
    # Cell goes CA1(n) --> CJ0(n))
    A[4N,2N] = k_i 
    # Cell goes CA2(n) --> CJ0(n))
    A[4N,3N] = k_i 

    return A
end

function ODEsolve_feedback_0mRNASTART(prms; 
    # Set maximum number of molecules
    N=150,
    timepoints = [0.0, 300.0, 500.0, 720.0, 1000, 2000, 4000],
    # Set time span
    tspan = (0.0,maximum(timepoints)),
    initialstate = :mRNAzero::Union{Symbol, Vector},
    )
    
    
    # function that creates a patrix from a set of parameters
    f = parameters -> A_matrix_feedback_0mRNASTART(parameters, N)

    # Function that multiplies the parameters
    F(x,p,t) = f(p)*x

    # Make initial state
    total_states = 4 
    
    if initialstate == :mRNAzero
        # All states are 0, but the CA0, which has 1 because all cells are in that state initially
        x0 = zeros(Float64, total_states*N+1);   
        # Now the initial state is the last state in the matrix
        x0[end] = 1.0;    
        x0 = x0./sum(x0)
    else 
        x0 = initialstate
    end
    
    
    masterProblem = ODEProblem(F,x0,tspan, prms)
    sol = solve(masterProblem, BS3(),saveat = timepoints);
    
    # You can use this bit of code to chech that the  
    # sum of the probabilities of boing in one of 
    # the individual stages is equal to 1 USEFUL FOR DEBUGGING !!!
     # probs = [sum(a[1:end]) for a in sol.u] ################
    
    # Calculate burst fraction
    
    burst_fract = [sum(a[N+1:3N]) for a in sol.u]
    
    mRNA = [a[1:N] .+ a[N+1:2N] .+ a[2N+1:3N] .+ a[3N+1:4N].+ append!([a[end]], zeros(N-1)) for a in sol.u]
    
    return mRNA, burst_fract, sol.t
    
end


