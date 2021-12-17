function _objective_function_sim(θ,X)

    # alias some parameters -
    α = θ[1]
    β = θ[2]
    R_firm = X[:,1]
    R_market = X[:,2]

    # compute the residual -
    ϵ = (R_firm .- α .- β*R_market)

    # return the total error -
    return sum(ϵ.^2)
end

function build_single_index_model(market::DataFrame, firm::DataFrame;
    risk_free_rate::Float64 = 5.15e-5)::SingleIndexModel


    # firm risk premimum -
    Y = (firm[!,:μ]) .- risk_free_rate

    # market risk premimum -
    M = (market[!,:μ]) .- risk_free_rate

    # compute the var of the market -
    β = cov(Y,M)/var(M)
    α = mean(Y) - (β)*mean(M)
    θ = [α,β]

	# compute the residuals -
    number_of_time_steps = length(Y)
    X = [ones(number_of_time_steps,1) M]
    r = X*θ - Y
    μ_R = mean(r)
    σ_R = std(r)

	# build model -
	model = SingleIndexModel()
	model.α = θ[1]
	model.β = θ[2]
	model.r = risk_free_rate
	model.ϵ = Normal(μ_R, σ_R)

    # return -
    return model
end

function build_single_index_model(market::DataFrame, firms::Dict{String,DataFrame}, 
    ticker_symbol_array::Array{String,1}, start::Date, stop::Date; risk_free_rate::Float64 = 5.15e-5)

    # initialize -
    model_dictionary = Dict{String, SingleIndexModel}()

    # extract out the date range from the market -
    market_df = extract_data_block_for_date_range(market, start, stop)

    # process -
    for ticker_symbol ∈ ticker_symbol_array

        # firm df -
        firm_df = extract_data_block_for_date_range(firms[ticker_symbol], start, stop);

        # data sets need to have the same number of rows -
        if (nrow(firm_df) == nrow(market_df))
            model = build_single_index_model(market_df, firm_df; risk_free_rate=risk_free_rate)
            model_dictionary[ticker_symbol] = model
        end
    end

    # return -
    return model_dictionary
end