function compute_contract_set_expiration(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # get some stuff the dictionary -
    contract_set = dictionary["contract_set"]
    simulation_price_dictionary = dictionary["simulated_price_range"]

    # setup asset price array -
    price_start = simulation_price_dictionary["start_price"]
    price_stop = simulation_price_dictionary["stop_price"]
    price_length = simulation_price_dictionary["length"]
	price_array = range(price_start,stop=price_stop,length=price_length) |> collect

    # load the contract set -
	r = build_simulation_contract_set(contract_set)
	if (isa(r.value,Exception) == true)
		return r.value
	end
    contract_set = r.value
    
    # simulate -
	r_pla = compute_option_profit_and_loss_at_expiration(contract_set,price_array)
	if (isa(r_pla.value,Exception) == true)
		return r_pla.value
	end
	pla = r_pla.value

    # return -
    return PSResult(pla)
end

function compute_put_price_binary_model(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # check: do we have contract set key?
    if (haskey(dictionary,"contract_set_parameters") == false)
        error = PSError("Missing the contract_set_parameters key in the request dictionary")
        return PSResult(error)
    end

    # get the contract set from the dictionary -
    contract_set = dictionary["contract_set_parameters"]
    r = build_simulation_contract_set(contract_set)
	if (isa(r.value,Exception) == true)
		return r
	end
    contract_dictionary = r.value
    ticker_symbol = contract_dictionary["ticker_symbol"]
    contract_sense_symbol = Symbol(contract_dictionary["sense"])
    expiration_date = Date(contract_dictionary["expiration"])
    contract_multiplier = contract_dictionary["contract_multiplier"]
    number_of_contracts = contract_dictionary["number_of_contracts"]
    premium_value = contract_dictionary["premium_value"]

    # do we have the binary simulation parameters?
    if (haskey(dictionary,"binary_simulation_parameters") == false)
        error = PSError("Missing the binary_simulation_parameters key in the request dictionary")
        return PSResult(error)
    end
    binary_lattice_dictionary = dictionary["binary_lattice_parameters"]
    
    # build binary lattice model -
    volatility = binary_lattice_dictionary["volatility"]
    riskFreeRate = binary_lattice_dictionary["riskFreeRate"]
    dividendRate = binary_lattice_dictionary["dividendRate"]
    numberOfLevels = binary_lattice_dictionary["numberOfLevels"]
    cal_days_to_exp = binary_lattice_dictionary["calendarDaysToExpiration"]
    days_to_expiration = (cal_days_to_exp)*(numberOfLevels/365)
	binary_lattice_model = PSBinaryLatticeModel(volatility, days_to_experiation,
        risk_free_rate, dividend_rate);

    # do we have an underlying price in the request?
    if (haskey(dictionary,"underlying_price_value") == false)
        error = PSError("Missing the underlying_price_value key in the request dictionary")
        return PSResult(error)
    end
    underlying_price_value = dictionary["underlying_price_value"]

    # get simulated strike price range -
    if (haskey(dictionary,"simulated_strike_price_dictionary") == false)
        error = PSError("Missing the simulated_strike_price_dictionary key in the request dictionary")
        return PSResult(error)
    end
    simulated_strike_price_dictionary = dictionary["simulated_strike_price_dictionary"]
    start_price = simulated_strike_price_dictionary["start_price"]
    stop_price = simulated_strike_price_dictionary["stop_price"]
    number_of_steps = simulated_strike_price_dictionary["number_of_steps"]

    # build the array -
    strike_price_array = collect(range(price_start,stop=price_stop,length=number_of_steps))

    # calc the contract cost for this set of strike prices
    predicted_price_array = Array{Float64,2}(undef,number_of_steps,4)
    for (index,strike_price) in enumerate(strike_price_array)

        # ok, first we create the contract object (in this case put, buy)
    	putOptionContract = PSPutOptionContract(ticker_symbol, expiration_date, 
            strike_price, 0.0, number_of_contracts; 
            sense=contract_sense_symbol, contractMultiplier=contract_multiplier)
    
        # create a contract set -
        contract_set = Set{PSAbstractAsset}()
        push!(contract_set, putOptionContract);

        # run the pricing calc -
        local result = option_contract_price(contract_set,binary_lattice_model,
            underlying_price_value; earlyExercise=true)
    
        if (isa(result.value,Exception) == true)
            return result
        end
        results_tuple = result.value;

        # grab the computed price -
        p = results_tuple.cost_calculation_result.option_contract_price_array[1];
        predicted_price_array[index,1] = strike_price
        predicted_price_array[index,2] = underlying_price_value
        predicted_price_array[index,3] = max(0,(strike_price - underlying_price_value))
        predicted_price_array[index,4] = p
    end

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