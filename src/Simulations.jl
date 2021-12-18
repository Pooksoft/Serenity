function simulate_sim_random_walk(model::SingleIndexModel, initial_price::Float64, 
    market_return::Array{Float64,1}; N::Int64 = 100, Δt = 1.0)

    # get the number of steps -
    number_of_steps = length(market_return)

    # initialize some storage -
    simulated_price_array = Array{Float64,2}(undef, number_of_steps, N)
    simulated_price_array[1,1:N] .= initial_price

    # get stuff from the model -
    α = model.α
    β = model.β
    r = model.r
    ϵ = model.ϵ

    # main loop -
    for path_index ∈ 1:N

        # what is the random pertrubation?
        ϵ_array = rand(ϵ, number_of_steps)
        
        # simulation loop -
        for simulation_step_index ∈ 2:number_of_steps
            
            # what is the return for the market at this time index?
            r_market = market_return[simulation_step_index - 1]

            # ok, compute the return in this time period -
            μ = α + β*(r_market - r) + ϵ_array[simulation_step_index] - r
        
            # now that we have the return μ, let's calc the new price -
            old_price = simulated_price_array[simulation_step_index - 1, path_index]
            new_price = old_price*(1+μ*Δt)
            simulated_price_array[simulation_step_index, path_index] = new_price
        end
    end

    # return -
    return simulated_price_array
end

# short cut methods -
(model::SingleIndexModel)(initial_price::Float64, market_return::Array{Float64,1}; Δt = 1.0) = 
    simulate_sim_random_walk(model, initial_price, market_return; Δt = Δt);