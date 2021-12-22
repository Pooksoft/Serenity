# ------------------------------------------------------------------------------------------------------------------------------- #
# compute_log_return_array(data::DataFrame,map::Pair{Symbol,Symbol}; Δt::Float64 = 1.0) -> DataFrame

# Computes μ of historical price data. 

# # Arguments:
# data_table  	DataFrame holding the historical price data
# map 	        Different data APIs return data with different field names. The map arg connects the time field to price field
# Δt 		    Time step size. Default is 1.0 in the units you want the μ to be calculated in
# ------------------------------------------------------------------------------------------------------------------------------ #
function compute_log_return_array(data_table::DataFrame, map::Pair{Symbol,Symbol}; Δt = (1.0 / 365.0))

    # initialize -
    (number_of_rows, _) = size(data_table)
    return_table = DataFrame(timestamp = Date[], P1 = Float64[], P2 = Float64[], μ = Float64[])

    # main loop -
    for row_index = 2:number_of_rows

        # grab the date -
        tmp_date = data_table[row_index, map.first]

        # grab the price data -
        yesterday_close_price = data_table[row_index-1, map.second]
        today_close_price = data_table[row_index, map.second]

        # compute the diff -
        μ = (1 / Δt) * log(today_close_price / yesterday_close_price)

        # push! -
        push!(return_table, (tmp_date, yesterday_close_price, today_close_price, μ))
    end

    # return -
    return return_table
end

function compute_log_return_array(ticker_symbol_array::Array{String,1}, data_tables::Dict{String,DataFrame}, 
    map::Pair{Symbol,Symbol}; Δt = (1.0 / 365.0))

    # initialize -
    number_of_ticker_symbols = length(ticker_symbol_array)
    data_dictionary = Dict{String,DataFrame}()

    # process -
    for ticker_symbol_index ∈ 1:number_of_ticker_symbols

        # get the ticker symbol for this index -
        ticker_symbol = ticker_symbol_array[ticker_symbol_index]

        # get the data frame for this index -
        df = data_tables[ticker_symbol]

        # compute the return -
        data_dictionary[ticker_symbol] = compute_log_return_array(df, map; Δt = Δt)
    end

    # return -
    return data_dictionary
end

function compute_fractional_return_array(data_table::DataFrame, map::Pair{Symbol,Symbol};
    multiplier::Float64=1.0)

    # initialize -
    (number_of_rows, _) = size(data_table)
    return_table = DataFrame(timestamp = Date[], P1 = Float64[], P2 = Float64[], μ = Float64[])

    # main loop -
    for row_index = 2:number_of_rows

        # grab the date -
        tmp_date = data_table[row_index, map.first]

        # grab the price data -
        yesterday_close_price = data_table[row_index-1, map.second]
        today_close_price = data_table[row_index, map.second]

        # compute the diff -
        μ = ((today_close_price/yesterday_close_price) - 1.0)*multiplier;

        # push! -
        push!(return_table, (tmp_date, yesterday_close_price, today_close_price, μ))
    end

    # return -
    return return_table
end

function compute_fractional_return_array(ticker_symbol_array::Array{String,1}, data_tables::Dict{String,DataFrame}, map::Pair{Symbol,Symbol}; 
    multiplier::Float64=1.0)

    # initialize -
    number_of_ticker_symbols = length(ticker_symbol_array)
    data_dictionary = Dict{String,DataFrame}()

    # process -
    for ticker_symbol_index ∈ 1:number_of_ticker_symbols

        # get the ticker symbol for this index -
        ticker_symbol = ticker_symbol_array[ticker_symbol_index]

        # get the data frame for this index -
        df = data_tables[ticker_symbol]

        # compute the return -
        data_dictionary[ticker_symbol] = compute_fractional_return_array(df, map; multiplier=multiplier)
    end

    # return -
    return data_dictionary
end

function compute_analytical_geometric_brownian_motion_trajectory(model::GeometricBrownianMotionModel, 
    initial_condition::Float64, number_of_steps::Int64; Δt::Float64 = 1.0, N::Int64 = 100)

    # initialize -
    state_array = Array{Float64,2}(undef, number_of_steps, N)
    state_array[1,:] .= initial_condition

    # get data from the model -
    μ = model.μ
    σ = model.σ

    # calculate the random term -
    d = Normal{Float64}(0.0, 1.0)

    # calculate the drift term -
    drift_term = μ*Δt

    # main: sample paths -
    for sample_path_index ∈ 1:N

        # generate new sequence of steps for this sample path -
        random_terms = σ*sqrt(Δt)*rand(d,number_of_steps)
        
        # solve the model -
        for time_step_index ∈ 2:number_of_steps
            
            old_state = state_array[time_step_index-1,sample_path_index]
            
            # compute the price difference -
            Δ = old_state*(drift_term+random_terms[time_step_index])

            # compute the new state -
            new_state = old_state+Δ

            # capture -
            state_array[time_step_index,sample_path_index] = new_state
        end
    end

    # return -
    return state_array
end

function compute_random_walk_model_trajectory(model::Distribution, initial_price::Float64,
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

function compute_rwm_cumulative_probabilty(compare::Function, price_array::Array{Float64,1})

    # initialize -
    number_of_samples = length(price_array)
    tmp_array = BitArray(undef, (number_of_samples, 1))

    # main -
    for sample_index = 1:number_of_samples

        # get the sample price -
        sample_price = price_array[sample_index]

        # check: which is larger, sample or target price?
        compare(sample_price) ? tmp_array[sample_index] = 1 : tmp_array[sample_index] = 0
    end

    # sum the tmp_array -
    number_of_larger_values = sum(tmp_array)

    # compute the probability -
    return (number_of_larger_values / number_of_samples)
end

function compute_rwm_cumulative_probabilty(price_array::Array{Float64,1}, target_price::Float64)

    # initialize -
    number_of_samples = length(price_array)
    tmp_array = Array{Int64,1}()

    # main -
    for sample_index = 1:number_of_samples

        # get the sample price -
        sample_price = price_array[sample_index]

        # check: which is larger, sample or target price?
        sample_price <= target_price ? push!(tmp_array, 1) : push!(tmp_array, 0)
    end

    # sum the tmp_array -
    number_of_larger_values = sum(tmp_array)

    # compute the probability -
    return (number_of_larger_values / number_of_samples)
end

function compute_minvar_portfolio_allocation(μ,Σ,target_return::Float64; 
    w_lower::Float64 = 0.0, w_upper::Float64 = 1.0)

    # initialize -
    number_of_assets = length(μ)
    w = Variable(number_of_assets)
    risk = quadform(w,Σ)
    ret  = dot(w,μ)

    # setup problem -
    p = minimize(risk)
    p.constraints += [sum(w)==1.0, w_lower <= w, w <= w_upper, ret >= target_return]
    Convex.solve!(p, SCS.Optimizer(verbose = false))

    # return -
    return (p.status, evaluate(w), p.optval, evaluate(ret))
end

function compute_cybernetic_portfolio_allocation(μ,Σ)

    # how many assets do we have?
    𝒫 = length(μ)

    # initialize -
	term_array = Array{Float64,1}(undef, 𝒫)
	for term_index ∈ 1:𝒫
		term_array[term_index] = max(0.0, μ[term_index]/Σ[term_index,term_index])
	end

    # compute the u-variable -
	𝒵 = sum(term_array)

    # if we have no good options (all the returns are negative, then all allocations would be zero)
    if (𝒵 == 0)
        u_variable_array = zeros(𝒫)
    else
        u_variable_array = (1/𝒵)*term_array
    end

    # return -
    return u_variable_array
end

function compute_average_fractional_return_and_covariance(tickers::Array{String,1}, data::Dict{String,DataFrame}, 
    start::Date, stop::Date)

    # need to grab the number of time steps -
    test_df = Serenity.extract_data_block_for_date_range(data[tickers[1]], start, stop)
	number_of_time_steps = nrow(test_df)

    # initialize -
    μ_bar = Array{Float64,1}()
    R_array = Array{Float64,2}(undef, number_of_time_steps, length(tickers))

    # process each ticker -
    for (ticker_index, ticker) ∈ enumerate(tickers)

        # get the data for the data range -
        df = data[ticker]
        df_slice = extract_data_block_for_date_range(df,start,stop)

        # compute the return -
		avg_val = mean(df_slice[!,:μ])
		push!(μ_bar, avg_val)

        # build R_array -
		for step_index = 1:number_of_time_steps
			R_array[step_index, ticker_index] = df_slice[step_index,:μ]
		end
    end

    # compute the covariance -
    Σ = cov(R_array)

    # return -
    return (μ_bar,Σ)
end


# short cut method RWMC method -
(model::Distribution)(initial_price::Float64, number_of_steps::Int64; number_of_sample_paths = 1) =
    compute_random_walk_model_trajectory(model, initial_price, number_of_steps;
        number_of_sample_paths = number_of_sample_paths)

