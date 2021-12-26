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

function simulate_insample_portfolio_allocation(tickers::Array{String,1}, 
    return_dictionary::Dict{String,DataFrame}, ω::Array{Float64,1}, initial_budget::Float64, 
    start::Date, stop::Date; multiplier::Float64 = 1.0)

    # initialize -
    df = extract_data_block_for_date_range(return_dictionary[tickers[1]], start, stop);
    number_of_timesteps = nrow(df)
    number_of_tickers = length(tickers)
    wealth_array = Array{Float64,2}(undef, number_of_timesteps, number_of_tickers)
    
    # initialize the wealth array -
    for ticker_index ∈ 1:number_of_tickers
        wealth_array[1,ticker_index] = ω[ticker_index]*initial_budget
    end

    for (index, ticker) ∈ enumerate(tickers)

        # grab -
        df = extract_data_block_for_date_range(return_dictionary[ticker], start, stop);

        for time_index ∈ 2:number_of_timesteps
            
            # get the return -
            r = multiplier*df[time_index, :μ];

            # compute the wealth -
            old_wealth = wealth_array[time_index-1,index]
            new_wealth = (1+r)*old_wealth
            wealth_array[time_index,index] = new_wealth
        end
    end

    # return the wealth_array -
    return wealth_array
end

function simulate_random_walk_model_trajectory(model::Distribution, initial_price::Float64,
    number_of_steps::Int64; number_of_sample_paths = 1)

    # initialize -
    number_of_steps = number_of_steps + 1
    price_array = Array{Float64,2}(undef, number_of_steps, number_of_sample_paths)

    # insert the first value -
    price_array[1, 1:number_of_sample_paths] .= initial_price

    # how many samples do we need?
    sample_return_array = rand(model, number_of_steps, number_of_sample_paths)

    # compute the price -
    for sample_path_index = 1:number_of_sample_paths
        for step_index = 2:number_of_steps
            price_array[step_index, sample_path_index] = price_array[step_index-1, sample_path_index] + sample_return_array[step_index, sample_path_index]
        end
    end

    # return -
    return price_array
end

# short cut methods -
(model::SingleIndexModel)(initial_price::Float64, market_return::Array{Float64,1}; Δt = 1.0) = 
    simulate_sim_random_walk(model, initial_price, market_return; Δt = Δt);

# short cut method RWMC method -
(model::Distribution)(initial_price::Float64, number_of_steps::Int64; number_of_sample_paths = 1) =
    simulate_random_walk_model_trajectory(model, initial_price, number_of_steps;
        number_of_sample_paths = number_of_sample_paths)
