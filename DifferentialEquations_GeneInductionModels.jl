module DifferentialEquations_GeneInductionModels

using DifferentialEquations
using ProgressMeter
using PyPlot
using Distributed


mutable struct ReactionCtes
    ka::Float64
    ki::Float64
    kon::Float64
    koff::Float64
    ktc::Float64
    kdeg::Float64
end

############################
## GENE EXPRESSION MODELS ##
############################
# 20th Jun 2018, v1.0 #
######################
"""
-----------------
| MODEL SPECIES |
-----------------
Cell states:
    - c_i0 Cell in a starting state where it cannot express Il12b
    - c_a  Cell in a state where it can express Il12b
    - c_i1 Cell in a second inactive state where it has been able to express Il12b, 
      but the negative feedback has been activated and it cannot longer express Il12b
Species inside the cell:
    - gen_off Number of Il12b alleles where transcription cannot occur
    - gen_on Number of Il12b allelles where transcription can occur
    - mRNA Number of messenger molecules of Il12b
"""
telegraph_model = @reaction_network TEL begin
        r_on, gen_off --> gen_on
        r_off, gen_on --> gen_off
        r_tc, gen_on --> gen_on + mRNA
        r_dg, mRNA --> 0
    end r_on r_off r_tc r_dg

induction_model = @reaction_network TEL_ind begin
        r_a, c_i0 --> c_a # Transition to active state
        0, c_a --> c_i0 # This should not happen
        r_on, c_a + gen_off --> c_a + gen_on
        r_off, gen_on --> gen_off
        r_tc, gen_on --> gen_on + mRNA
        r_dg, mRNA --> 0
    end r_a r_on r_off r_tc r_dg

feedback_model = @reaction_network TEL_ind_feed begin
        r_a, c_i0 --> c_a # Transition to active state
        0, c_a --> c_i0 # This should not happen
        r_i, c_a --> c_i1 # Negative feedback is active in the cell
        0, c_i1 --> c_a # This neither
        r_on, c_a + gen_off --> c_a + gen_on
        r_off, gen_on --> gen_off
        r_tc, gen_on --> gen_on + mRNA
        r_dg, mRNA --> 0
    end r_a r_i r_on r_off r_tc r_dg


###################
## MODEL SOLVERS ##
###################

function telegraph_solve(k::ReactionCtes; ploidity = 2, time_span = (0.0,7200.0), 
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
    start_cell = [ploidity, 0, mRNA]
    model = telegraph_model
    react_ctes = (k.kon, k.koff, k.ktc, k.kdeg)
    
    #Define Discrete problems (We are dealing with absolute number of molecules)
    prob = DiscreteProblem(start_cell, time_span, react_ctes)
    jump_prob = JumpProblem(prob,Direct(),model)
    # Solve
    # The saveat might be useful to seed up the simulations
    if saveat != 0
        return sol = solve(jump_prob, FunctionMap(), saveat = saveat, save_everystep = save_everystep)
    else
        return sol = solve(jump_prob, FunctionMap())
    end
end


function induction_solve(k::ReactionCtes; ploidity = 2, time_span = (0.0,7200.0), 
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
    model = induction_model
    react_ctes = (k.ka, k.kon, k.koff, k.ktc, k.kdeg)
    
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


function feedback_solve(k::ReactionCtes; ploidity = 2, time_span = (0.0,7200.0), 
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
    start_cell = [1, 0, 0, ploidity, 0, mRNA]
     model = feedback_model
     react_ctes = (k.ka, k.ki, k.kon, k.koff, k.ktc, k.kdeg)
    
    #Define Discrete problems (We are dealing with absolute number of molecules)
     prob = DiscreteProblem(start_cell, time_span, react_ctes)
     jump_prob = JumpProblem(prob,Direct(),model)
    # Solve
    # The saveat might be useful to seed up the simulations
    if saveat != 0
        return sol = solve(jump_prob, FunctionMap(), saveat = saveat,  save_everystep = save_everystep, dense = false)
    else
        return sol = solve(jump_prob, FunctionMap())
    end
end

####################################################################
## FUNCTIONS TO GET SIMULATION AND REAL DATA FOR DISTANCE MEASURE ##
####################################################################


function get_simulation_vector(cellctes, model_solve; 
        saveat = 1800, 
        n_simulations = 300, 
        maxtime = 7200,
        timepoints = [0.0, 3600.0, 5400.0, 7200.0],
        mRNA = true,
        save_everystep = false,
        BurstFraction = true
    )
     simcell = model_solve(cellctes, saveat = saveat,  save_everystep = save_everystep)
        

    tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
    times = simcell.t[tim_ind]
    mRNAs = Array{Int64, 2}(length(tim_ind), n_simulations)
    # Quantify the number of mRNA molecules present in the cell
    mRNAs[:, 1] = [x[end] for x in simcell.u[tim_ind]]
    burst_cells =  Array{Int64, 2}(length(tim_ind), n_simulations) 
    
   
    # Check wether the cell is bursting at those timepoints
     v_length = length(simcell.u[tim_ind][1])
     burst_cells[:, 1] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    
    for cell in 2:n_simulations
        simcell = model_solve(cellctes, saveat = saveat,  save_everystep = save_everystep)
        tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
        times = simcell.t[tim_ind]
        mRNAs[:, cell] = [x[end] for x in simcell.u[tim_ind]]
        burst_cells[:, cell] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    end
            
            
    burst_fract = [sum(burst_cells[i, :])/size(burst_cells)[2] for i in 1:size(burst_cells)[1]]
    
    if mRNA == true && BurstFraction == true
       return [mRNAs[i, :] for i in 1:size(mRNAs)[1]], burst_fract
    end
    
end


function get_simulation_vector(cellctes1, cellctes2, model_solve; 
        saveat = 1800, 
        n_simulations = 300, 
        maxtime = 7200,
        timepoints = [0.0, 3600.0, 5400.0, 7200.0],
        mRNA = true,
        save_everystep = false,
        BurstFraction = true
    )
    simcell = model_solve(cellctes1, cellctes2, saveat = saveat,  save_everystep = save_everystep)
    tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
    times = simcell.t[tim_ind]
    mRNAs = Array{Int64, 2}(length(tim_ind), n_simulations)
    # Quantify the number of mRNA molecules present in the cell
    mRNAs[:, 1] = [x[end] for x in simcell.u[tim_ind]]
    burst_cells =  Array{Int64, 2}(length(tim_ind), n_simulations) 
    
   
    # Check wether the cell is bursting at those timepoints
    v_length = length(simcell.u[tim_ind][1])
    burst_cells[:, 1] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    
    for cell in 2:n_simulations
        simcell = model_solve(cellctes1, cellctes2, saveat = saveat,  save_everystep = save_everystep)
        tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
        times = simcell.t[tim_ind]
        mRNAs[:, cell] = [x[end] for x in simcell.u[tim_ind]]
        burst_cells[:, cell] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    end
            
            
    burst_fract = [sum(burst_cells[i, :])/size(burst_cells)[2] for i in 1:size(burst_cells)[1]]
    
    if mRNA == true && BurstFraction == true
       return [mRNAs[i, :] for i in 1:size(mRNAs)[1]], burst_fract
    end
    
end



function get_simulation_vector_p(cellctes, model_solve; 
        saveat = 1800, 
        n_simulations = 300, 
        maxtime = 7200,
        timepoints = [0.0, 3600.0, 5400.0, 7200.0],
        mRNA = true,
        BurstFraction = true
    )
    simcell = model_solve(cellctes, saveat = saveat)
    tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
    times = simcell.t[tim_ind]
    # They will be shared by the different processes so I need to create a shared array
    mRNAs = SharedArray{Int64, 2}(length(tim_ind), n_simulations)
    # Quantify the number of mRNA molecules present in the cell
    mRNAs[:, 1] = [x[end] for x in simcell.u[tim_ind]]
    # They will be shared by the different processes so I need to create a shared array
    burst_cells =  SharedArray{Int64, 2}(length(tim_ind), n_simulations) 
    
   
    # Check wether the cell is bursting at those timepoints
    v_length = length(simcell.u[tim_ind][1])
    burst_cells[:, 1] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    
     @sync @distributed for cell in 2:n_simulations
        simcell = model_solve(cellctes, saveat = saveat)
        tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
        times = simcell.t[tim_ind]
        mRNAs[:, cell] = [x[end] for x in simcell.u[tim_ind]]
        burst_cells[:, cell] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    end
            
            
    burst_fract = [sum(burst_cells[i, :])/size(burst_cells)[2] for i in 1:size(burst_cells)[1]]
    
    if mRNA == true && BurstFraction == true
       return [mRNAs[i, :] for i in 1:size(mRNAs)[1]], burst_fract
    end
    
end
        
function get_simulation_vector_pmap(cellctes, model_solve; 
        saveat = 1800, 
        n_simulations = 300, 
        maxtime = 7200,
        timepoints = [0.0, 3600.0, 5400.0, 7200.0],
        mRNA = true,
        BurstFraction = true
    )
            
    simcell = model_solve(cellctes, saveat = saveat)
    tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
    times = simcell.t[tim_ind]
    # They will be shared by the different processes so I need to create a shared array
    mRNAs = SharedArray{Int64, 2}(length(tim_ind), n_simulations)
    # Quantify the number of mRNA molecules present in the cell
    mRNAs[:, 1] = [x[end] for x in simcell.u[tim_ind]]
    # They will be shared by the different processes so I need to create a shared array
    burst_cells =  SharedArray{Int64, 2}(length(tim_ind), n_simulations) 
    
   
    # Check wether the cell is bursting at those timepoints
    v_length = length(simcell.u[tim_ind][1])
    burst_cells[:, 1] =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
    
    function simcell(nsim)
        simcell = model_solve(cellctes, saveat = saveat)
        tim_ind = [something(findfirst(isequal(ii), simcell.t), 0) for ii in timepoints]
        times = simcell.t[tim_ind]
        mRNA = [x[end] for x in simcell.u[tim_ind]]
        burst_cell =  [x[v_length - 1] >= 1 for x in simcell.u[tim_ind]]
        return mRNA, burst_cell
    end
                    
     res = pmap(simcell, 1:(n_simulations-1))
     mRNA_res = [i[1] for i in res]
     burst_res = [i[2] for i in res]
                    
    for i in 1:length(res)
        mRNAs[:, i+1] = mRNA_res[i]
        burst_cells[:, i+1] = burst_res[i]
    end
            
    
    burst_fract = [sum(burst_cells[i, :])/size(burst_cells)[2] for i in 1:size(burst_cells)[1]]
    
    if mRNA == true && BurstFraction == true
       return [mRNAs[i, :] for i in 1:size(mRNAs)[1]], burst_fract
    end
    
end

function get_simulation_vector_from_realdata(exp_data, genotype; 
        saveat = 1800, 
        n_simulations = 100, 
        maxtime = 7200,
        timepoints = [0, 60, 90, 120],
        mRNA = true,
        BurstFraction = true          
    )
            
    f(x) = [i == genotype for i in x[:Genotype]]
    table_data = exp_data[f(exp_data), :]
    
    mRNA_v = Array{Array, 1}(length(timepoints))
    burst_v = Array{Float64, 1}(length(timepoints))
    
    i = 0
    for tim in timepoints
        i +=1
        f(x, time) = [i == time for i in x[:Timepoint]]
        mRNA_v[i] = table_data[f(table_data, tim), :N_exon]
        burst_f = [i>=1 for i in table_data[f(table_data, tim), :total_TS_Cell]]
        burst_v[i] =  sum(burst_f)/length(burst_f)
    end
            
    if mRNA == true && BurstFraction == true
       return [i for i in mRNA_v], burst_v
    end
         
            
end
        
############################################
## helper functions for distance measures ##
############################################
function norm_frequency(df::Array)
    # It returns an array with the normalized frequencies from 0 to the maximum in the array
    samples_in = length(df)
    maximum_in = Int(maximum(df))
    new_array = zeros(maximum_in + 1)
    
    for i in 0:maximum_in
        
        new_array[i+1] = count(x->x==i,df)/samples_in
    end
  
    return new_array
end

function cumnorm_frequency(df::Array)
    # It returns an array with the cumulative normalized frequencies from 0 to the maximum in the array
    df = norm_frequency(df)
    for i in 2:length(df)
        df[i] = df[i] + df[i-1]
    end
    df
end

function KSstatistic(pop1, pop2)
    pop1 = cumnorm_frequency(pop1)
    pop2 = cumnorm_frequency(pop2)
    l1 = length(pop1)
    l2 = length(pop2)
        
    if l2 > l1
        append!(pop1, fill(1.0, l2-l1))
    elseif l1 > l2
        append!(pop2, fill(1.0, l1-l2))
    end
        
    #ks = maximum(abs(pop2 - pop1))  # This worked previously
    ks = maximum(abs.(pop2 - pop1))  
end
                        
function gini(wages, wagefrequencies)
    Swages_previous = wages[1]*wagefrequencies[1]
    Gwages = Swages_previous*wagefrequencies[1]
    @inbounds for i = 2:length(wages)
        freq = wagefrequencies[i]
        Swages_current = Swages_previous + wages[i]*freq
        Gwages += freq * (Swages_current+Swages_previous)
        Swages_previous = Swages_current
    end
    return 1.0 - Gwages/Swages_previous
end

function gini(d::Array{Int64,1})
    freqs = norm_frequency(d)
    counts = collect(0:1:maximum(d))
    g = gini(counts, freqs)
    if isnan(g)
        return 1
    else
        return g
    end
end

function gini_abs_diff(
        sim1::Tuple{Array{Array{Int64,1},1},Array{Float64,1}}, 
        sim2::Tuple{Array{Array{Int64,1},1},Array{Float64,1}}
    )
    g = 0
    # Do not take 0 timepoint into consideration
    for i in 2:length(sim1[1])
        g = g + abs(gini(sim1[1][i]) - gini(sim2[1][i]))
    end
    return g
end

#######################
## DISTANCE MEASURES ##
#######################

function KSstatistics_plus_squared_diff_cellsbursting(
        data1::Tuple{Array{Array{Int64,1},1},Array{Float64,1}}, 
        data2::Tuple{Array{Array{Int64,1},1},Array{Float64,1}})
    # Get the number of mRNA per cell
    mRNA1 = data1[1]
    mRNA2 = data2[1]
    # Calculate the KS statistic
    KSs = [KSstatistic(mRNA1[i], mRNA2[i]) for i in 1:length(mRNA1)]
    # Get the fraction of cells with a TS
    b1 = data1[2]
    b2 = data2[2]
    # Calculate the square difference
    squared_diff = [abs(b1[i]-b2[i])^2 for i in 1:length(b1)]
    
    distm = sum(KSs) + sum(squared_diff)
end
    
function KSstatistics_plus_diff_cellsbursting(
        data1::Tuple{Array{Array{Int64,1},1},Array{Float64,1}}, 
        data2::Tuple{Array{Array{Int64,1},1},Array{Float64,1}})
    # Get the number of mRNA per cell
    mRNA1 = data1[1]
    mRNA2 = data2[1]
    # Calculate the KS statistic
    KSs = [KSstatistic(mRNA1[i], mRNA2[i]) for i in 1:length(mRNA1)]
    # Get the fraction of cells with a TS
    b1 = data1[2]
    b2 = data2[2]
    # Calculate the square difference
    squared_diff = [abs(b1[i]-b2[i]) for i in 1:length(b1)]
    
    distm = sum(KSs) + sum(squared_diff)
end

function KSstatistics_plus_diff_cellsbursting_gini(
        data1::Tuple{Array{Array{Int64,1},1},Array{Float64,1}}, 
        data2::Tuple{Array{Array{Int64,1},1},Array{Float64,1}})
    # Get the number of mRNA per cell
    mRNA1 = data1[1]
    mRNA2 = data2[1]
    # Calculate the KS statistic
    KSs = [KSstatistic(mRNA1[i], mRNA2[i]) for i in 2:length(mRNA1)]
    # Get the fraction of cells with a TS
    b1 = data1[2]
    b2 = data2[2]
    # Calculate the square difference
    squared_diff = [abs(b1[i]-b2[i]) for i in 2:length(b1)]
    # Calculate gini diferences
    ginidiff = gini_abs_diff(data1, data2)
    # Calculate distance
    distm = sum(KSs) + sum(squared_diff) + ginidiff
end



########################
## PLOTTING FUNCTIONS ##
########################

function plot_activealleles(virtualcell; model_solve = telegraph_solve, color = "blue", cellname = "Cell", alpha = 1)
    if model_solve == telegraph_solve
        data = [x[2] for x in virtualcell.u]
    elseif model_solve == induction_solve
        data = [x[4] for x in virtualcell.u]
    elseif model_solve == feedback_solve
        data = [x[5] for x in virtualcell.u]
    end
    t = virtualcell.t
    PyPlot.plot(t, data, label = cellname, color = color, lw = 2, alpha = alpha)
    ylabel("Active alleles")
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false) # Hide the right edge of the axis
    axes[:xaxis][:set_ticks_position]("bottom")
    axes[:yaxis][:set_ticks_position]("left")
    ylim(0, 2.5)
    xticks([], [])
    yticks([ 1, 2])
end

function plot_cellstate(virtualcell; model_solve = telegraph_solve, color = "blue", cellname = "Cell")
    if model_solve == telegraph_solve
        data = [1 for x in virtualcell.u]
    elseif model_solve == induction_solve
        data = [if x[1] == 1 0 else 1 end for x in virtualcell.u]
    elseif model_solve == feedback_solve
        data = [if x[1] == 1 0 elseif x[2] == 1 1 else 2 end  for x in virtualcell.u]
    end
    t = virtualcell.t
    PyPlot.plot(t, data, label = cellname, color = color, lw = 2)
    ylabel("Cell State")
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false) # Hide the right edge of the axis
    axes[:xaxis][:set_ticks_position]("bottom")
    axes[:yaxis][:set_ticks_position]("left")
    ylim(-0.3, 2.5)
    xticks([], [])
    yticks([ 0, 1, 2], ["I0", "A", "I1"])
end

function plot_mRNA(virtualcell; model_solve = telegraph_solve, color = "blue", cellname = "Cell", alpha = 1)
    data = [x[end] for x in virtualcell.u]
    t = virtualcell.t
    PyPlot.plot(t, data, label = cellname, color = color, lw = 2, alpha = alpha)
    ylabel("mRNA")
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false) # Hide the right edge of the axis
    axes[:xaxis][:set_ticks_position]("bottom")
    axes[:yaxis][:set_ticks_position]("left")
    xticks([], [])
end

################################################################
## PLOTTING VIRTUAL CELL POPUTLATION SCREENSHOT DISTRIBUTIONS ##
################################################################
    
function mRNA_distributions_overtime(cellctes, model_solve; saveat = 10, n_simulations = 1000, maxtime = 7200)
    simcell = model_solve(cellctes, saveat = saveat)
    tim_ind = [findfirst( simcell.t, i) for i in 0:saveat:maximum(simcell.t)]
    times = simcell.t[tim_ind]
    mRNAs = Array{Int64, 2}(length(tim_ind), n_simulations)
    mRNAs[:, 1] = [x[end] for x in simcell.u[tim_ind]]

    p = Progress(n_simulations-1, 1)
    
    for cell in 2:n_simulations
        simcell = model_solve(cellctes, saveat = saveat)
        tim_ind = [findfirst( simcell.t, i) for i in 0:saveat:maximum(simcell.t)]
        times = simcell.t[tim_ind]
        mRNAs[:, cell] = [x[end] for x in simcell.u[tim_ind]]
        next!(p)
    end
        
    mRNAs, times
    
end

function plot_distributions(dist_time; columns = true, plot_times = collect(0.0:1000:7200.0), 
        xlim_ = maximum(dist_time[1]),
        ylim_ = size(dist_time[1])[2]/3, 
        color = "blue"
    )
    ind = [findfirst(dist_time[2], i) for i in plot_times]
    
    cnt = 0
    for i in ind
        cnt += 1
        subplot(1, length(ind), cnt) 
        d = dist_time[1][i, :]
        plt[:hist](d, color = color)
        t = dist_time[2][i]
        title("time = $t")
        xlim(-1, xlim_)
        ylim(0, ylim_)

    
    end
    
   plt[:tight_layout]()
end

end
