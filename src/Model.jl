# For the moment, let's hardcode the CRR movement function - not sure how we can encode a "callback"
function _movement(latticeModel::PSBinaryLatticeModel)
    
    # compute -
    volatility = latticeModel.volatility
    timeToExercise = latticeModel.timeToExercise
    riskFreeRate = latticeModel.riskFreeRate
    numberOfLevels = latticeModel.numberOfLevels

    # what is the dT?
    Δt = (timeToExercise/numberOfLevels)

    # compute u -
    U = exp(volatility * √Δt)
    D = 1 / U

    # return -
    return (U,D)
end

function compute_equity_price_binary_model(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # initialize -  
    # ...  

    # check: do we have binary lattice parameters in the dictionary?
    if (haskey(dictionary,"binary_lattice_parameters") == false)
        error = PSError("Missing the binary_lattice_parameters key in the request dictionary")
        return PSResult(error)
    end
    binary_lattice_dictionary = dictionary["binary_lattice_parameters"]
    volatility = binary_lattice_dictionary["volatility"]
    risk_free_rate = binary_lattice_dictionary["riskFreeRate"]
    dividend_rate = binary_lattice_dictionary["dividendRate"]
    numberOfLevels = binary_lattice_dictionary["numberOfLevels"]
    cal_days_to_exp = binary_lattice_dictionary["calendarDaysToExpiration"]
    days_to_expiration = (cal_days_to_exp)*(1/365)
	binary_lattice_model = PSBinaryLatticeModel(volatility, days_to_expiration,
        risk_free_rate, dividend_rate);

    # check: do we have an initial underlying price value?
    if (haskey(dictionary,"underlying_price_value") == false)
        error = PSError("Missing the underlying_price_value key in the request dictionary")
        return PSResult(error)
    end
    underlying_price_value = dictionary["underlying_price_value"]

    # compute prob distribution for this underlying 0
    latticeModel = binary_lattice_model
    movementFunction = _movement
    baseUnderlyingPrice = underlying_price_value
    binomial_price_array = Array{Float64,2}(undef,(cal_days_to_exp+1),3)
    ti_array = range(0,stop=cal_days_to_exp,step=1) |> collect
    for (index,timeStepIndex) in enumerate(ti_array)

        # call method on options kit -
        result = compute_underlying_price_distribution(timeStepIndex, latticeModel, baseUnderlyingPrice, movementFunction)
        if (isa(result.value,Exception) == true)
            return result
        end
        pd = result.value

        # compute the expected price -
        expected_price = sum(pd[:,3].*pd[:,4])

        # compute the std for the binomial price -
        term_1 = sum(pd[:,3].*(pd[:,4].^2))
        std_binomial_price = sqrt(term_1 .- expected_price.^2)

        # capture -
        binomial_price_array[index,1] = index
        binomial_price_array[index,2] = expected_price
        binomial_price_array[index,3] = std_binomial_price
    end

    # return -
    return PSResult(binomial_price_array)
end