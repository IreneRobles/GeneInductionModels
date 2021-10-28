using Distributions
using DifferentialEquations
try import LinearAlgebra end


"""
Function to evaluate the log likelihood of the data given a particular set of parameters
"""

function loglikelihood(arrayofprobabilities, countdata)
    arrayofprobabilities = [if a == 0 eps(a) else a end for a in arrayofprobabilities]
    sum([log(abs.(arrayofprobabilities[a+1])) for a in countdata]) 
end


# Functions to run the MCMC metropolis algorithm

function mcmc_metropolis(
        countdata,
        timepoints,
        initial_guess, # Intitial guess of what the parameters can be
        model_solve, 
        jump_proposal, # Size of the step that modifies the parameters each time
        Lchain; # Number of jumps
        
        # Options
        burn::Integer = 500,
        step::Integer = 500,
        printFreq::Integer = 100,
        priors::Union{Array, Symbol} = :none, # It does not work
        initialstate = :mRNAzero,
        acceptanceratio = 1,
        fix_degradation = false,
        update_priors = true,
        update_proposal = true,
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
    else
        priors = [Truncated(Normal(a, sqrt(a)), 0 , Inf) for a in old_guess]
        loglik_mRNA = sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])])
        
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
        else
            newloglik_mRNA = 0.0
            for (ip, prr) in enumerate(priors)
                  newloglik_mRNA +=log(pdf(prr, new_guess[ip]))
            end
            if !isinf(newloglik_mRNA); 
                 newloglik_mRNA += sum([loglikelihood(ode_[1][a], countdata[1][a]) for a in 1:length(ode_[1])]); 
            end
            
            
        end
        
        likelihoodlog[ii] = loglik_mRNA
        # Likelihood ratio between the two points
        ratio = newloglik_mRNA - loglik_mRNA 
        
        if  ratio > -rand(acceptanceratio)[1]  # Accept the point with probability min(a,1)
            old_guess = new_guess
            # This is something I made up
            if update_priors == true
                priors = [Truncated(Normal(a, a*2), 0 , Inf) for a in old_guess]
            end
            # This is something I made up
            if update_proposal == true
                jump_proposal = MvNormal(0.01.*sqrt.(old_guess))
            end
            loglik_mRNA = newloglik_mRNA
            acc +=1
        end
        
        chain[ii, :] = old_guess
        
        if printFreq != 0 
            if ii % printFreq == 0
                println("Completed iteration $ii out of $Lchain, $acc particles accepted \n")
            end
        end

    end
    
    accrat = round(acc/Lchain, 2)
    println("Acceptance ratio of $accrat (%acc out of %Lchain)")
    chainRed = chain[burn:step:end, :]
    liks = likelihoodlog[burn:step:end] 
    
    return chainRed, liks
    
    
end