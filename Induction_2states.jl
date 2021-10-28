induction_model2 = @reaction_network TEL_ind begin
        r_a, c_i --> c_a # Transition to active state
        0, c_a --> c_i # This should not happen
    # Reactions for inactive state
        r_on, c_i + gen_off --> c_i + gen_on
        r_off, c_i + gen_on --> gen_off + c_i
        r_tc, c_i + gen_on --> gen_on + mRNA + c_i
        r_dg, c_i + mRNA --> c_i
    # Reactions for active state
        r_on2, c_a + gen_off --> c_a + gen_on
        r_off2, c_a + gen_on --> gen_off + c_a
        r_tc2, c_a + gen_on --> gen_on + mRNA + c_a
        r_dg2, c_a + mRNA --> c_a
    
end r_a r_on r_off r_tc r_dg r_on2 r_off2 r_tc2 r_dg2


function induction_solve2(
        k1::DifferentialEquations_GeneInductionModels.ReactionCtes, 
        k2::DifferentialEquations_GeneInductionModels.ReactionCtes; 
        ploidity = 2, time_span = (0.0,7200.0), 
    init_mRNA = -1, # If mRNA is not change, cell will start in an approximation of the steady state
    saveat = 0,
    save_everystep = true
    )
    # Get the mean number of mRNA molecules
    if init_mRNA < 0
        mRNA = 0
    else 
        mRNA = init_mRNA
    end
    
    
    # This vector represents the inital state of the cell
    start_cell = [1, 0, ploidity, 0, mRNA]
    model = induction_model2
    react_ctes = (k1.ka, k1.kon, k1.koff, k1.ktc, k1.kdeg, k2.kon, k2.koff, k2.ktc, k2.kdeg)
    
    #Define Discrete problems (We are dealing with absolute number of molecules)
    prob = DiscreteProblem(start_cell, time_span, react_ctes)
    jump_prob = JumpProblem(prob,Direct(),model)
    # Solve
    # The saveat might be useful to seed up the simulations
    if saveat != 0
        return sol = solve(jump_prob, FunctionMap(), saveat = saveat,   save_everystep = save_everystep)
    else
        return sol = solve(jump_prob, FunctionMap())
    end

end

function A_matrix_induction2(parameters, N)
    
    k_a = parameters[1]
    
    k1_on = parameters[2]
    k1_off = parameters[3]
    k1_tc = parameters[4]
    k1_dg = parameters[5]
    
    k2_on = parameters[6]
    k2_off = parameters[7]
    k2_tc = parameters[8]
    k2_dg = parameters[9]
    
    
    # 3 possible states 0, 1, or , 2 alleles activated
    total_states = 6 
    
    A = zeros(Float64, total_states*N, total_states*N)

    ## Assemble A row by row

    #####################################################################################

    for ii = 1:N
        
        ### For CI0 #########################################
        # ii - 1 because the ii correspond to 0 molecules
        # Cell goes out of CI0
        A[ii,ii] = -(2*k1_on+k1_dg*(ii-1) + k_a)
        # Cell goes CI0(n+1) --> CI0(n)
        if ii != N
            A[ii,ii+1] = k1_dg*ii #
        end
        # Cell goes CI1(n) --> CI0(n)
        A[ii,ii+N] = k1_off #
        #####################################################
        
        ### For CI1 #########################################
        # Cell goes CI0(n) --> CI1(n)
        A[ii+N,ii] = 2*k1_on #
        # Cell goes out of CI1
        if ii != N
            A[ii+N,ii+N] = -(k1_tc+k1_off+k1_on+k1_dg*(ii-1)+ k_a)
        else
            A[ii+N,ii+N] = -(k1_off+k1_on+k1_dg*(ii-1)+ k_a)
        end
        # Cell goes CI1(n+1) --> CI1(n)
        if ii != N
            A[ii+N,ii+N+1] = k1_dg*ii # 
        end
        # Cell goes CI1(n-1) --> CI1(n)
        if ii != 1
            A[ii+N,ii+N-1] = k1_tc #
        end
        # Cell goes CI2(n) --> CI1(n)
        A[ii+N, ii+2N] = 2*k1_off #
        #####################################################
        
        ### For CI2 #########################################
        # Cell goes CI1(n) --> CI2(n)
        A[ii+2N,ii+N] = k1_on #
        # Cell goes out of CI2
        if ii != N
            A[ii+2N,ii+2N] = -(2*k1_tc+2*k1_off+k1_dg*(ii-1)+ k_a) #
        else
            A[ii+2N,ii+2N] = -(2*k1_off+k1_dg*(ii-1)+ k_a)
        end
        # Cell goes CA2(n+1) --> CA2(n)
        if ii != N
            A[ii+2N,ii+2N+1] = k1_dg*ii #
        end
        # Cell goes CA2(n-1) --> CA2(n)
        if ii != 1
            A[ii+2N,ii+2N-1] = 2*k1_tc #
        end
        
        #####################################################
        #####################################################
        #####################################################
        
        
        ### For CA0 #########################################
        # ii - 1 because the ii correspond to 0 molecules
        # Cell goes out of CA0
        A[ii+3N,ii+3N] = -(2*k2_on+k2_dg*(ii-1))
        # Cell goes CI0 --> CA0
        A[ii+3N,ii] = k_a
        # Cell goes CA0(n+1) --> CA0(n)
        if ii != N
            A[ii+3N,ii+1+3N] = k2_dg*ii #
        end
        # Cell goes CA1(n) --> CA0(n)
        A[ii+3N,ii+4N] = k2_off #
        #####################################################
        
        ### For CA1 #########################################
        # Cell goes CA0(n) --> CA1(n)
        A[ii+4N,ii+3N] = 2*k2_on #
        # Cell goes CI1 --> CA1
        A[ii+4N,ii+N] = k_a
        # Cell goes out of CA1
        if ii != N
            A[ii+4N,ii+4N] = -(k2_tc+k2_off+k2_on+k2_dg*(ii-1))
        else
            A[ii+4N,ii+4N] = -(k2_off+k2_on+k2_dg*(ii-1))
        end
        # Cell goes CA1(n+1) --> CA1(n)
        if ii != N
            A[ii+4N,ii+4N+1] = k2_dg*ii # 
        end
        # Cell goes CA1(n-1) --> CA1(n)
        if ii != 1
            A[ii+4N,ii+4N-1] = k2_tc #
        end
        # Cell goes CA2(n) --> CA1(n)
        A[ii+4N, ii+5N] = 2*k2_off #
        #####################################################
        
        ### For CA2 #########################################
        # Cell goes CA1(n) --> CA2(n)
        A[ii+5N,ii+4N] = k2_on #
        # Cell goes CI2 --> CA2
        A[ii+5N,ii+2N] = k_a
        # Cell goes CA2(n) --> CA2(n+ 1) or CA2(n) --> CA1(n) or CA2(n) --> CA2(n - 1)
        if ii != N
            A[ii+5N,ii+5N] = -(2*k2_tc+2*k2_off+k2_dg*(ii-1)) #
        else
            A[ii+5N,ii+5N] = -(2*k2_off+k2_dg*(ii-1))
        end
        # Cell goes CA2(n+1) --> CA2(n)
        if ii != N
            A[ii+5N,ii+5N+1] = k2_dg*ii #
        end
        # Cell goes CA2(n-1) --> CA2(n)
        if ii != 1
            A[ii+5N,ii+5N-1] = 2*k2_tc #
        end
        #####################################################
        
    end
    
    return A
end

function ODEsolve_induction2(prms; 
    # Set maximum number of molecules
    N=10,
    timepoints = [0.0, 3000.0, 5000.0, 7200.0, 10000, 20000, 40000],
    # Set time span
    tspan = (0.0,maximum(timepoints)),
    initialstate = :mRNAzero::Union{Symbol, Vector},
    )
    
    
    # function that creates a patrix from a set of parameters
    f = parameters -> A_matrix_induction2(parameters, N)

    # Function that multiplies the parameters
    F(x,p,t) = f(p)*x

    # Make initial state
    total_states = 6 
    
    if initialstate == :mRNAzero
        # All states are 0, but the CA0, which has 1 because all cells are in that state initially
        x0 = zeros(Float64, total_states*N);   
        # Now the initial state is the last state in the matrix
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
    # return  probs = [sum(a[1:end]) for a in sol.u] ################
    
    # Calculate burst fraction
    
    burst_fract = [sum(a[N+1:3N]) + sum(a[4N+1:6N]) for a in sol.u]
    
   

    
    mRNA = [a[1:N] .+ a[N+1:2N] .+ a[2N+1:3N] .+ a[3N+1:4N] .+ a[4N+1:5N] .+ a[5N+1:6N] for a in sol.u]
    
    return mRNA, burst_fract, sol.t
    
end