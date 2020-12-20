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