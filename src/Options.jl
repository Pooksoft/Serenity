function compute_option_contract_set_expiration(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # get some stuff the dictionary -
    # do we have an underlying price in the request?
    if (haskey(dictionary,"contract_set_parameters") == false)
        error = PSError("Missing the contract_set_parameters key in the request dictionary")
        return PSResult(error)
    end

    # load the contract set -
	r = build_simulation_contract_set(dictionary)
	if (isa(r.value,Exception) == true)
		return r.value
	end
    contract_set = r.value

     # do we have an underlying price in the request?
     if (haskey(dictionary,"underlying_price_value") == false)
        error = PSError("Missing the underlying_price_value key in the request dictionary")
        return PSResult(error)
    end
    underlying_price_value = dictionary["underlying_price_value"]
    price_array = Array{Float64,1}()
    push!(price_array, underlying_price_value)
    
    # simulate -
	r_pla = compute_option_profit_and_loss_at_expiration(contract_set, price_array)
	if (isa(r_pla.value,Exception) == true)
		return r_pla.value
	end
    pla = r_pla.value
    
    # setup return -
    simulated_pl_dict = Dict{String,Any}()
    simulated_pl_dict["underlying_price_value"] = underlying_price_value
    simulated_pl_dict["profit_loss_value"] = pla[2]

    # return -
    return PSResult(simulated_pl_dict)
end

function compute_put_price_binary_model(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # check: do we have contract set key?
    if (haskey(dictionary,"contract_set_parameters") == false)
        error = PSError("Missing the contract_set_parameters key in the request dictionary")
        return PSResult(error)
    end

    # get the contract set from the dictionary -
    contract_dictionary = dictionary["contract_set_parameters"][1]
    ticker_symbol = contract_dictionary["ticker_symbol"]
    contract_sense_symbol = Symbol(contract_dictionary["sense"])
    expiration_date = Date(contract_dictionary["expiration"])
    contract_multiplier = contract_dictionary["contract_multiplier"]
    number_of_contracts = contract_dictionary["number_of_contracts"]
    contract_strike_price_value = contract_dictionary["contract_strike_price"]
    premium_value = 0.0

    # do we have the binary simulation parameters?
    if (haskey(dictionary,"binary_lattice_parameters") == false)
        error = PSError("Missing the binary_lattice_parameters key in the request dictionary")
        return PSResult(error)
    end
    binary_lattice_dictionary = dictionary["binary_lattice_parameters"]
    
    # build binary lattice model -
    volatility = binary_lattice_dictionary["volatility"]
    risk_free_rate = binary_lattice_dictionary["riskFreeRate"]
    dividend_rate = binary_lattice_dictionary["dividendRate"]
    numberOfLevels = binary_lattice_dictionary["numberOfLevels"]
    cal_days_to_exp = binary_lattice_dictionary["calendarDaysToExpiration"]
    days_to_expiration = (cal_days_to_exp)*(1/365)
	binary_lattice_model = PSBinaryLatticeModel(volatility, days_to_expiration,
        risk_free_rate, dividend_rate);

    # do we have an underlying price in the request?
    if (haskey(dictionary,"underlying_price_value") == false)
        error = PSError("Missing the underlying_price_value key in the request dictionary")
        return PSResult(error)
    end
    underlying_price_value = dictionary["underlying_price_value"]

    # ok, first we create the contract object (in this case put, buy)
    putOptionContract = PSPutOptionContract(ticker_symbol, expiration_date, 
        contract_strike_price_value, premium_value, number_of_contracts; 
        sense=contract_sense_symbol, contractMultiplier=contract_multiplier)

    # create a contract set -
    contract_set = Set{PSAbstractAsset}()
    push!(contract_set, putOptionContract);

    # run the pricing calc -
    result = option_contract_price(contract_set, binary_lattice_model,
        underlying_price_value; earlyExercise=true)
    if (isa(result.value, Exception) == true)
        return result
    end
    results_tuple = result.value

    # get the price -
    p = results_tuple.cost_calculation_result.option_contract_price_array[1];
    contract_value = p
    if (p === nothing)
        contract_value = 0.0
    end
    
    # we are going to include some additional stuff in the 
    predicted_price_dict = Dict{String,Any}()
    predicted_price_dict["contract_strike_price"] = contract_strike_price_value
    predicted_price_dict["underlying_price_value"] = underlying_price_value
    predicted_price_dict["intrinsic_contract_value"] = max(0,(contract_strike_price_value - underlying_price_value))
    predicted_price_dict["total_contract_value"] = contract_value   

    # return -
    return PSResult(predicted_price_dict)
end

function compute_contract_price_binary_model(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # initialize -
    
    # check: do we have contract set key?
    if (haskey(dictionary,"contract_set") == false)
        error = PSError("Missing the contract_set key in the request dictionary")
        return PSResult(error)
    end

    # get the contract set from the dictionary -
    contract_set = dictionary["contract_set"]
    r = build_simulation_contract_set(contract_set)
	if (isa(r.value,Exception) == true)
		return r.value
	end
    contract_set = r.value

    # do we have an underlying price in the request?
    if (haskey(dictionary,"underlying_price_value") == false)
        error = PSError("Missing the underlying_price_value key in the request dictionary")
        return PSResult(error)
    end
    underlying_price_value = dictionary["underlying_price_value"]

    # do we have the binary simulation parameters?
    if (haskey(dictionary,"binary_simulation_parameters") == false)
        error = PSError("Missing the binary_simulation_parameters key in the request dictionary")
        return PSResult(error)
    end

    # get the binary simulation parameters dictionary, and associated parameters -
    binary_sim_parameters_dictionary = dictionary["binary_simulation_parameters"]

    @show binary_sim_parameters_dictionary

    volatility = binary_sim_parameters_dictionary["volatility"]
    timeToExercise = binary_sim_parameters_dictionary["timeToExercise"]
    riskFreeRate = binary_sim_parameters_dictionary["drift_rate"]
    dividendRate = binary_sim_parameters_dictionary["dividendRate"]
    numberOfLevels = binary_sim_parameters_dictionary["numberOfLevels"]
    lattice_model = PSBinaryLatticeModel(volatility,timeToExercise,riskFreeRate,dividendRate; 
        numberOfLevels=numberOfLevels)

    # call the options calculation -
    compute_result = option_contract_price(contract_set, lattice_model,underlying_price_value;
        earlyExercise=true,numberOfLevels=numberOfLevels)

    # wrap and return -
    return PSResult(compute_result)
end