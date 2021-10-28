using PyPlot

function compare_gillespie_ode_burst(gillespie, ode_result; c_gillespie = "purple", c_ode = "orange")
    timepoints = ode_result[3]
    
    burst1 = gillespie[2]
    burst2 = ode[2]
    
    scatter(timepoints, burst1, c = c_gillespie, label = "Gillespie")
    plot(timepoints, burst1, c = c_gillespie)
    
    scatter(timepoints, burst2, c = c_ode, label = "ODE")
    plot(timepoints, burst2, c = c_ode)
    
    ylabel("Fraction of cells with at least 1 gene ON")
    
    xlabel("Simulation time")
    
    ax = gca()
    ax[:spines]["right"][:set_visible](false)
    ax[:spines]["top"][:set_visible](false)
    
    legend()
    title("Cells with at least 1 gene ON")
end

function compare_gillespie_ode_mRNAs(gillespie, ode_result; c_gillespie = "purple", c_ode = "orange")
    
    figure(figsize = (5*(length(ode_result[3])-1), 5))
    
    mRNA1 = [DifferentialEquations_GeneInductionModels.norm_frequency(a) for a in gillespie[1]]
    mRNA2 = ode[1]
    

    for a  in 2:length(ode_result[3])
        
    subplot(1, length(ode_result[3])-1, a-1)
    
     t = ode_result[3][a]
    
    scatter(0:length(mRNA1[a])-1, mRNA1[a], c = c_gillespie, label = "Gillespie")
    plot(0:length(mRNA1[a])-1, mRNA1[a], c = c_gillespie , linewidth = 4)
    
    scatter(0:length(mRNA2[a])-1, mRNA2[a], c = c_ode, label = "ODE")
    plot(0:length(mRNA2[a])-1, mRNA2[a], c = c_ode)
    
    ylabel("Fraction of cells")
    
    xlabel("mRNA molecules")
    xlim(0 - 0.5, length(mRNA1[a]) + 0.5)
    
    ax = gca()
    ax[:spines]["right"][:set_visible](false)
    ax[:spines]["top"][:set_visible](false)
    
    
    title("t = $t")
        
    end
    
    legend()
    plt[:tight_layout]()
end
