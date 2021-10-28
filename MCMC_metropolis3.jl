using Distributions
using DifferentialEquations
# This fix problems related to julia versions 0.6 and 0.7
try import LinearAlgebra catch end


# This function convert vectors into reaction constants, useful for test_mcmc
function turn_params_in_ReactionCtes(vector, ode_model)
    if ode_model == ODEsolve_telegraph_model
        return DifferentialEquations_GeneInductionModels.ReactionCtes(0, 0, vector...)
        elseif ode_model == ODEsolve_telegraph_steady_state
            return DifferentialEquations_GeneInductionModels.ReactionCtes(0, 0, vector...)
        elseif  ode_model == ODEsolve_induction_0mRNASTART
            return DifferentialEquations_GeneInductionModels.ReactionCtes(vector[1], 0, vector[2:end]...)
        elseif  ode_model == ODEsolve_feedback_0mRNASTART_withJ1J2
            return DifferentialEquations_GeneInductionModels.ReactionCtes(vector...)
    end
    
end


"""
Function to evaluate the log likelihood of the data given a particular set of parameters
"""
function loglikelihood(arrayofprobabilities, countdata)
    # Fix issues with 0
    arrayofprobabilities = [if a == 0 eps() else a end for a in arrayofprobabilities]
    sum([log(abs.(arrayofprobabilities[a+1])) for a in countdata]) 
end

function loglikelihood_TSS(ode, counts)
    # get the number of cells in each timepoint
    ncells = [length(ii) for ii in counts[1]]
    # Find probability in ODE to be not bursting and bursting, store in a vector
    ode_states  = [[1-ii, ii] for ii in  ode[2]] 
    # Fix issues with 0
    for state in 1:length(ode_states)
        for a in 1:length(ode_states[state])
            if ode_states[state][a] == 0 
                ode_states[state][a] = eps()
            end 
        end
    end
    # Find out how many cells in the data are bursting and not bursting
    cells_states = round.(counts[2].*ncells, 0)
    cells_states = [[ncells[ii] - cells_states[ii], cells_states[ii]] for ii in 1:length(ncells)]
    
    return sum([log(ode_states[ii][1])*cells_states[ii][1] + log(ode_states[ii][2])*cells_states[ii][2] for ii in 1:length(ncells)])
    

end


# Functions to run the MCMC metropolis algorithm

function mcmc_metropolis(
        countdata,
        timepoints,
        initial_guess, # Intitial guess of what the parameters can be
        model_solve, # ode function generator
        jump_proposal, # Size of the step that modifies the parameters each time
        Lchain; # Number of jumps
        
        # Options
        burn::Integer = 500,
        step::Integer = 500,
        printFreq::Integer = 100,
        priors::Union{Array, Symbol} = :none, # It does not work
        initialstate = :mRNAzero,
        acceptanceratio = 1, # Keep this as 1, other is a bad idea
        fix_degradation = false,
        update_priors = false,
        update_proposal = true,
        proposal_stdfactor = 0.01
    )
    
    ## Make ProgressBar ## For impatience people like me
    #p = Progress(Lchain, dt = 1)
    
    old_guess = initial_guess
    n_param = length(old_guess)
    chain = zeros(Float64, Lchain, n_param)
    likelihoodlog = zeros(Float64, Lchain)
    chain[1, :] = old_guess
    acc = 0
    
    Nmax = maximum([maximum(a) for a in countdata[1]])+1
    
    ode_ =  model_solve(old_guess, 
        N=Nmax, 
        timepoints = timepoints, 
        initialstate = initialstate)
    
    if priors == :none
        loglik_mRNA = sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])])
        loglik_mRNA += loglikelihood_TSS(ode_, countdata)
    else
        if priors == :standar
            priors = [Truncated(Normal(a, sqrt(a)*2), 0 , Inf) for a in old_guess]
        else 
            priors = priors
        end
        loglik_mRNA = sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])])
         loglik_mRNA += loglikelihood_TSS(ode_, countdata)
        for (ip, prr) in enumerate(priors)
               loglik_mRNA +=log(pdf(prr, old_guess[ip]))
        end
    end
    likelihoodlog[1] = loglik_mRNA
    
    for ii = 2:Lchain
        #next!(p)
        # Obtain new guess from proposal distribution
        new_guess = abs.(old_guess + rand(jump_proposal)) # Make sure that this numbers are positive
        if fix_degradation == true
            new_guess[end] = old_guess[end] # Last parameter (kdeg is fix)
        end
        ode_ =  model_solve(new_guess, N=Nmax, 
                    timepoints = timepoints, 
                    initialstate = initialstate)
        if priors == :none
             newloglik_mRNA = sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])]);      
            newloglik_mRNA += loglikelihood_TSS(ode_, countdata)
        else
            newloglik_mRNA = 0.0
            for (ip, prr) in enumerate(priors)
                  newloglik_mRNA +=log(pdf(prr, new_guess[ip]))
            end
            if !isinf(newloglik_mRNA); 
                newloglik_mRNA += sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])]); 
                newloglik_mRNA += loglikelihood_TSS(ode_, countdata)
            end
            
            
        end
        
        # Likelihood ratio between the two points
        ratio = newloglik_mRNA - loglik_mRNA 
        
        if  Base.MathConstants.e^ratio > rand(1)[1]*acceptanceratio  # Accept the particel with a certain probability 
            old_guess = new_guess
            """
            # This is a bad idea according to Rowan (and he is right)
            if update_priors == true
                priors = [Truncated(Normal(a, a*2), 0 , Inf) for a in old_guess]
            end
            """
            if update_proposal == true
                jump_proposal = MvNormal(proposal_stdfactor.*sqrt.(old_guess))
            end
            loglik_mRNA = newloglik_mRNA
            acc +=1
        end
        
        likelihoodlog[ii] = loglik_mRNA
        
        chain[ii, :] = old_guess
        
        if printFreq != 0 
            if ii % printFreq == 0
                println("Completed iteration $ii out of $Lchain, $acc particles accepted \n")
            end
        end
        
       
        
            
    end
    
    accrat = round(acc/Lchain, digits=4)
    println("Acceptance ratio of $accrat (%acc out of %Lchain) \n")
    chainRed = chain[burn:step:end, :]
    liks = likelihoodlog[burn:step:end] 
    return chainRed, liks
    
    
end

function test_mcmc(; 
        gillespie_simulation::Union{String, Function} = "",
        ode_model = "",
        timepoints = [0.0, 60.0, 90.0, 120.0],
        TestName = "",
        gillespie_parameters = "",
        initial_guess = "",
        n_cells = 1000,
        Lchain = 1000,
        jump_proposal = [],
        acceptanceratio = "",
        printFreq = 100,
        burn = 1000,
        fix_degradation = false,
        update_priors = true,
        update_proposal = true,
        proposal_stdfactor = 0.01
    
    )
    
    # Save MCMC information
    
    info = Dict()
    info[:Name] = TestName
    info[:ODEModel] = "$ode_model"
    info[:Timepoints] = timepoints
    info[:GillespieSimulation] = "$gillespie_simulation"
    info[:TrueParameters] = gillespie_parameters
    info[:InitialGuess] = initial_guess
    info[:GillespieCells] = n_cells
    
    
    

    filename = "MCMC__" *TestName *"__"* string(Dates.now()) *"__MCMC_OUT.jld"
    
    # Simulate smFISH data
    
    
    cellctes = turn_params_in_ReactionCtes(gillespie_parameters, ode_model)
    if length(timepoints) > 1
        saveat = minimum(abs.([timepoints[v] - timepoints[v-1] for v in 2:length(timepoints)]))
    else 
        saveat = timepoints[1]
    end
    t_sim = @time begin
        println("Gillespie simulation time:")
        sim = DifferentialEquations_GeneInductionModels.get_simulation_vector(cellctes, gillespie_simulation; 
            saveat = saveat,    
            n_simulations = n_cells, 
                timepoints = timepoints,
                mRNA = true,
                BurstFraction = true
            )
    end
    
    info[:GillespieTime] = t_sim
    
    Nmax = maximum([maximum(a) for a in sim[1]])+1
    
    info[:Nmax] = Nmax
    
    # Save timings (One never knows when they will be useful)
    
    t_ode1 = @time begin
    
        println("ODE simulation time 1:")

        # Perform ODE simulation with true parameters

        ode = ode_model(gillespie_parameters; 
            N=Nmax,
            initialstate = :mRNAzero,
            timepoints = timepoints
            )
    end;
    
    info[:ODETime_1] = t_ode1
    
    t_ode2 = @time begin
    
        println("ODE simulation time 2:")

        # Perform ODE simulation with true parameters

        ode = ode_model(gillespie_parameters; 
            N=Nmax,
            initialstate = :mRNAzero,
            timepoints = timepoints
            )
    end;
    
    info[:ODETime_2] = t_ode2
    
    
    info[:Lchain] = Lchain
    info[:AcceptanceRatio] = acceptanceratio
    info[:fix_degradation] = fix_degradation
    info[:update_priors] = update_priors
    info[:update_proposal] = update_proposal
    info[:burn] = burn
    info[:step]= 500
    info[:proposal_stdfactor] = proposal_stdfactor
    
    jump_proposal = MvNormal(proposal_stdfactor.*sqrt.(initial_guess))
    priors = [Truncated(Normal(a, a), 0 , Inf) for a in initial_guess]
    
    t_mcmc = @time begin
    
    
    mcmc_out = mcmc_metropolis(
        sim,
        timepoints,
        initial_guess, # Intitial guess of what the parameters can be
        ode_model, 
        jump_proposal, # Size of the step that modifies the parameters each time
        Lchain; # Number of jumps
        
        # Options
        burn = burn,
        step= 500,
        printFreq= printFreq,
        priors = priors,
        initialstate = :mRNAzero,
        acceptanceratio = acceptanceratio,
        fix_degradation = fix_degradation,
        update_priors = update_priors,
        update_proposal = update_proposal,
        proposal_stdfactor = proposal_stdfactor
    )
    
    end
    
    info[:MCMC] = t_ode2
    
    JLD.save(filename, "INFO", info, "PARAMS", mcmc_out[1], "LOGLIKELIHOOD", mcmc_out[2])
    
    
    
    
end
