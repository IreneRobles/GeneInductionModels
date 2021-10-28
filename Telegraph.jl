# Function that generates the transition probability matrix
# Very prout to say that today, 24 Oct 2018, the matrix seem to work

"""
# This version is longer in code lines because it codes separately the states with 0, N molecules and all the other
# So on 29th October I decided to make a shorter version after making the first models 
# I did it because otherwise it would have been a hell to debug


function A_matrix_telegraph(parameters, N)
    
    k_on = parameters[1]
    k_off = parameters[2]
    k_tc = parameters[3]
    k_dg = parameters[4]
    
    
    # 3 possible states 0, 1, or , 2 alleles activated
    total_states = 3 
    
    A = zeros(Float64, total_states*N, total_states*N)

    ## Assemble A row by row
    
    #####################################################################################
    ### Special case is when the cell has 0 molecules of mRNA (it cannot degrade any) ###
    #####################################################################################
    ### For CA0
    # Cell goes CA0(0) --> CA1(0) 2k_on since cell is doploid
    A[1,1] = -2*k_on #
    # Cell goes CA0(1) --> CA0(n) 
    A[1,2] = k_dg #
    # Cell goes CA1(0) --> CA0(0)
    A[1,1+N] = k_off #
    
    ### For CA1
    # Cell goes CA0(n) --> CA1(n)
    A[1+N,1] = 2*k_on #
    # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n)
    A[1+N,1+N] = -(k_off + k_tc + k_on) #
    # Cell goes CA1(n + 1) --> CA1(n)
    A[1+N,2+N] = k_dg #
    # Cell goes CA2(n) --> CA1(n)
    A[1+N,1+2N] = 2*k_off #
    
    ### For CA2
    # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
    A[1+2N,1+2N] = -(2*k_tc + 2*k_off)
    # Cell goes CA2(n + 1) --> CA2(n)
    A[1+2N,2+2N] = k_dg #
    # Cell goes CA1(n) --> CA2(n)
    A[1+2N,1+N] = k_on #
    #####################################################################################

    for ii = 2:N-1
        ### For CA0
        # Cell goes CA0(n) --> CA1(n) or CA0(n) --> CA0(n - 1) 
        # ii - 1 because the ii correspond to 0 molecules
        A[ii,ii] = -(2*k_on+k_dg*(ii-1))
        # Cell goes CA0(n+1) --> CA0(n)
        A[ii,ii+1] = k_dg*ii #
        # Cell goes CA1(n) --> CA0(n)
        A[ii,ii+N] = k_off #
        
        ### For CA1
        # Cell goes CA0(n) --> CA1(n)
        A[ii+N,ii] = 2*k_on #
        # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n) or CA1(n) --> CA1(n - 1)
        A[ii+N,ii+N] = -(k_tc+k_off+k_on+k_dg*(ii-1))
        # Cell goes CA1(n+1) --> CA1(n)
        A[ii+N,ii+N+1] = k_dg*ii # 
        # Cell goes CA1(n-1) --> CA1(n)
        A[ii+N,ii+N-1] = k_tc #
        # Cell goes CA2(n) --> CA1(n)
        A[ii+N, ii+2N] = 2*k_off #
        
        
        ### For CA2
        # Cell goes CA1(n) --> CA2(n)
        A[ii+2N,ii+N] = k_on #
        # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
        A[ii+2N,ii+2N] = -(2*k_tc+2*k_off+k_dg*(ii-1)) #
        # Cell goes CA2(n+1) --> CA2(n)
        A[ii+2N,ii+2N+1] = k_dg*ii #
        # Cell goes CA2(n-1) --> CA2(n)
        A[ii+2N,ii+2N-1] = 2*k_tc #
        
    end
    
    
    #####################################################################################
    ### Special case is when the cell has maximum molecules of mRNA (it cannot degrade any) ###
    #####################################################################################
   
    ### For CA0
    # Cell goes CA0(n) --> CA1(n) or CA0(n) --> CA0(n - 1) 
    A[N,N] = -(2*k_on+k_dg*(N-1))
    # Cell goes CA1(n) --> CA0(n)
    A[N,2*N] = k_off
    
            
    ### For CA1
    # Cell goes CA0(n) --> CA1(n)
    A[2*N,N] = 2*k_on
    # Cell goes CA1(n) --> CA1(n+ 1) or CA1(n) --> CA0(n) or CA1(n) --> CA2(n) or CA1(n) --> CA1(n - 1)
    A[2*N,2*N] = -(k_off+k_on+k_dg*(N-1))
    # Cell goes CA1(n-1) --> CA1(n)
    A[2*N,2*N-1] = k_tc
    # Cell goes CA2(n) --> CA1(n)
    A[2*N, 3*N] = 2*k_off
    
    ### For CA2
    # Cell goes CA1(n) --> CA2(n)
    A[3*N,2*N] = k_on
    # Cell goes  CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
    A[3*N,3*N] = -(2*k_off+k_dg*(N-1))
    # Cell goes CA2(n-1) --> CA2(n)
    A[3*N,3*N-1] = 2*k_tc

    return A
end
"""
function A_matrix_telegraph(parameters, N)
    
    k_on = parameters[1]
    k_off = parameters[2]
    k_tc = parameters[3]
    k_dg = parameters[4]
    
    
    # 3 possible states 0, 1, or , 2 alleles activated
    total_states = 3 
    
    A = zeros(Float64, total_states*N, total_states*N)

    ## Assemble A row by row

    #####################################################################################

    for ii = 1:N
        ### For CA0 #########################################
        # ii - 1 because the ii correspond to 0 molecules
        # Cell goes out of CA0
        A[ii,ii] = -(2*k_on+k_dg*(ii-1))
        # Cell goes CA0(n+1) --> CA0(n)
        if ii != N
            A[ii,ii+1] = k_dg*ii #
        end
        # Cell goes CA1(n) --> CA0(n)
        A[ii,ii+N] = k_off #
        #####################################################
        
        ### For CA1 #########################################
        # Cell goes CA0(n) --> CA1(n)
        A[ii+N,ii] = 2*k_on #
        # Cell goes out of CA1
        if ii != N
            A[ii+N,ii+N] = -(k_tc+k_off+k_on+k_dg*(ii-1))
        else
            A[ii+N,ii+N] = -(k_off+k_on+k_dg*(ii-1))
        end
        # Cell goes CA1(n+1) --> CA1(n)
        if ii != N
            A[ii+N,ii+N+1] = k_dg*ii # 
        end
        # Cell goes CA1(n-1) --> CA1(n)
        if ii != 1
            A[ii+N,ii+N-1] = k_tc #
        end
        # Cell goes CA2(n) --> CA1(n)
        A[ii+N, ii+2N] = 2*k_off #
        #####################################################
        
        ### For CA2 #########################################
        # Cell goes CA1(n) --> CA2(n)
        A[ii+2N,ii+N] = k_on #
        # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
        if ii != N
            A[ii+2N,ii+2N] = -(2*k_tc+2*k_off+k_dg*(ii-1)) #
        else
            A[ii+2N,ii+2N] = -(2*k_off+k_dg*(ii-1))
        end
        # Cell goes CA2(n+1) --> CA2(n)
        if ii != N
            A[ii+2N,ii+2N+1] = k_dg*ii #
        end
        # Cell goes CA2(n-1) --> CA2(n)
        if ii != 1
            A[ii+2N,ii+2N-1] = 2*k_tc #
        end
        #####################################################
        
    end
    
    return A
end



function ODEsolve_telegraph_steady_state(prms; 
    # Set maximum number of molecules
    N=150,
    timepoints = [0.0, 300.0, 500.0, 720.0, 1000, 2000, 4000],
    # Set time span
    tspan = (0.0,maximum(timepoints)),
    initialstate = :mRNAzero::Union{Symbol, Vector},
    )
    # FInd Steady State
    Sol = LinearAlgebra.nullspace(A_matrix_telegraph(prms, N))
    # Divide by sum(Sol) or the possibility of being in each state will be more than 1
    Sol = Sol ./ sum(Sol)
    
    burst_fract = [sum(Sol[N+1:end])]
    mRNA = [Sol[1:N] .+ Sol[N+1:2N] .+ Sol[2N+1:3N]]
    
    return mRNA, burst_fract
end



function ODEsolve_telegraph_model(prms; 
    # Set maximum number of molecules
    N=150,
    timepoints = [0.0, 300.0, 500.0, 720.0, 1000, 2000, 4000],
    # Set time span
    tspan = (0.0,maximum(timepoints)),
    initialstate = :mRNAzero::Union{Symbol, Vector},
    )
    
    
    # function that creates a patrix from a set of parameters
    f = parameters -> A_matrix_telegraph(parameters, N)

    # Function that multiplies the parameters
    F(x,p,t) = f(p)*x

    # Make initial state
    total_states = 3 
    
    if initialstate == :mRNAzero
        # All states are 0, but the CA0, which has 1 because all cells are in that state initially
        x0 = zeros(Float64, total_states*N);    
        x0[1] = 1.0;    
        x0 = x0./sum(x0)
    else 
        x0 = initialstate
    end
    
    masterProblem = ODEProblem(F,x0,tspan, prms)
    sol = solve(masterProblem, BS3(),saveat = timepoints);
    
    # You can use this bit of code to chech that the  
    # sum of the probabilities of boing in one of 
    # the individual stages is equal to 1 USEFUL FOR DEBUGGING !!!
    probs = [sum(a[1:end]) for a in sol.u] ################
    
    # Calculate burst fraction
    
    burst_fract = [sum(a[N+1:end]) for a in sol.u]
    
    mRNA = [a[1:N] .+ a[N+1:2N] .+ a[2N+1:3N] for a in sol.u]
    
    return mRNA, burst_fract, sol.t
    
end