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

    # check: do we have the same number of rows in the two df?
    if (nrow(market) != nrow(firm))
        return nothing
    end 

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

function build_single_index_model(market::DataFrame, firm::DataFrame, 
    start::Date, stop::Date; risk_free_rate::Float64 = 5.15e-5)::SingleIndexModel

    # extract out the date range from the market -
    market_df = extract_data_block_for_date_range(market, start, stop);
    
    # extract out the firm df -
    firm_df = extract_data_block_for_date_range(firm, start, stop);

    # check: if we don't have the same number of rows, then return nothing
    if (nrow(market_df) != nrow(firm_df))
        return nothing
    end

    # return -
    return build_single_index_model(market_df, firm_df; risk_free_rate = risk_free_rate);
end


function build_single_index_model(market::DataFrame, firms::Dict{String,DataFrame}, 
    ticker_symbol_array::Array{String,1}, start::Date, stop::Date; 
        risk_free_rate::Float64 = 5.15e-5)

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

function build_geometric_brownian_motion_model(asset_return_array::Array{Float64,1})::GeometricBrownianMotionModel

    # initialize -
    gbm = GeometricBrownianMotionModel()
    
    # we need to estimate the mean and std from the return array -
    gbm.μ = mean(asset_return_array)
    gbm.σ = std(asset_return_array)

    # build the model -
    return gbm
end

function add_parameters_to_url_query_string(base::String, options::Dict{String,Any})::String

    # init -
    url_string = base

    parameters = ""
    for (key, value) in options
        parameters *= "$(key)=$(value)&"
    end

    # cut off trailing &
    query_parameters = parameters[1:end-1]

    # return -
    return url_string * query_parameters
end

function build(base::String, model::PolygonAggregatesEndpointModel; 
    apiversion::Int64 = 2)::String

    # get data from the API call data -
    adjusted = model.adjusted
    limit = model.limit
    sortdirection = model.sortdirection
    apikey = model.apikey
    ticker = model.ticker
    to = model.to
    from = model.from
    multiplier = model.multiplier
    timespan = model.timespan

    # build up the base string -
    base_url = "$(base)/v$(apiversion)/aggs/ticker/$(ticker)/range/$(multiplier)/$(timespan)/$(from)/$(to)?"

    # what keys are passed as parameters?
    options_dictionary = Dict{String,Any}()
	options_dictionary["adjusted"] = adjusted
	options_dictionary["sort"] = sortdirection
	options_dictionary["limit"] = limit
	options_dictionary["apiKey"] = apikey

    # return -
    return add_parameters_to_url_query_string(base_url, options_dictionary)
end

function build(base::String, model::PolygonOptionsContractReferenceEndpoint; 
    apiversion::Int64 = 3)::String

    # what are the fields for this api model?
    api_model_fieldnames = fieldnames(PolygonOptionsContractReferenceEndpoint);

    # build -
    base_url = "$(base)/v$(apiversion)/reference/options/contracts?"

    # build an options dictionary -
    options_dictionary = Dict{String,Any}()
    for fieldname ∈ api_model_fieldnames
        
        value = getproperty(model,fieldname)
        if (isnothing(value) == false)
            options_dictionary[String(fieldname)] = value
        end
    end

    # add some defaults to the options dictionary -
    get!(options_dictionary,"order","asc")
    get!(options_dictionary,"sort","strike_price")
    get!(options_dictionary,"limit", 1000) # default to max if not specified

    # now - add the parameters ...
    return add_parameters_to_url_query_string(base_url, options_dictionary)
end

function build_option_ticker_symbol(underlying::String, expiration::Date, type::PolygonOptionContractType, 
    K::Float64)

    # compute the ticker string -
    
    # get all the components -
    ticker_component = uppercase(underlying)
    YY = year(expiration) - 2000 # hack to get a two digit year 
    MM = lpad(month(expiration),2,"0")
    DD = lpad(day(expiration),2,"0")

    # what is the option contract type?
    contract_type = "C"
    if (type == put)
        contract_type = "P"
    end

    # compute the price code -
    strike_component = lpad(convert(Int64, K*1000), 8, "0");

    # build the string -
    ticker_string = "$(ticker_component)$(YY)$(MM)$(DD)$(contract_type)$(strike_component)"

    # return the ticker string -
    return ticker_string;
end