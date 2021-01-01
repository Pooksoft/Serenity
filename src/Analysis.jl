function compute_linear_return_equity(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # initialize -  
    # ...  

    try 
    
        # get the equity_price_array from the dictionary -
        if (haskey(dictionary,"equity_price_array") == false)
            error = PSError("Missing the equity_price_array key in the request dictionary")
            return PSResult(error)
        end
        equity_price_array_tmp = dictionary["equity_price_array"]
        equity_price_array = map(x->convert(Float64,x),equity_price_array_tmp)

        # run the computation -
        result = compute_linear_return_array(equity_price_array)
        if (isa(result.value, Exception) == true)
            return result
        end
        return_array = result.value

        # return -
        return PSResult(return_array)

    catch error
        return PSResult(error)
    end
end

function compute_log_return_equity(dictionary::Dict{String,Any})::PooksoftBase.PSResult

    # initialize -  
    # ...  

    try
        
        # get the equity_price_array from the dictionary -
        if (haskey(dictionary,"equity_price_array") == false)
            error = PSError("Missing the equity_price_array key in the request dictionary")
            return PSResult(error)
        end
        equity_price_array_tmp = dictionary["equity_price_array"]
        equity_price_array = map(x->convert(Float64,x),equity_price_array_tmp)
        
        # get log base - 
        if (haskey(dictionary,"base") == false)
            error = PSError("Missing the base key in the request dictionary")
            return PSResult(error)
        end
        log_base = dictionary["base"] # for now: we ignore this parameter, but we'll update and refactor

        # run the computation -
        result = compute_log_return_array(equity_price_array)
        if (isa(result.value, Exception) == true)
            return result
        end
        return_array = result.value

        # return -
        return PSResult(return_array)
    catch error
        return PSResult(error)
    end
end