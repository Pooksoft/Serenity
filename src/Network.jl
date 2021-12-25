function http_get_call_with_url(url::String)::Some

    try

        # should we check if this string is formatted as a URL?
        if (occursin("https://", url) == false)
            throw(ArgumentError("url $(url) is not properly formatted"))
        end

        # ok, so we are going to make a HTTP GET call with the URL that was passed in -
        response = HTTP.request("GET", url)

        # ok, so let's check if we are getting a 200 back -
        if (response.status == 200)
            return Some(String(response.body))
        else
            # create an error, and throw it back to the caller -
            throw(ErrorException("http status flag $(response.status) was returned from url $(url)"))
        end
    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function polygon_error_handler(request_body_dictionary::Dict{String, Any})::Tuple

    # initialize -
    error_response_dictionary = Dict{String,Any}()
    
    # what are my error keys?
    error_keys = [
        "status", "error", "request_id"
    ]
    for key ∈ error_keys
        error_response_dictionary[key] = request_body_dictionary[key]
    end

    # return -
    return (error_response_dictionary, nothing)
end

function process_aggregates_polygon_call_response(body::String)

    # convert to JSON -
    request_body_dictionary = JSON.parse(body)

    # before we do anything - check: do we have an error?
    status_flag = request_body_dictionary["status"]
    if (status_flag == "ERROR")
        return polygon_error_handler(request_body_dictionary)
    end

    # initialize -
    header_dictionary = Dict{String,Any}()
    df = DataFrame(

        volume=Float64[],
        volume_weighted_average_price=Float64[],
        open=Float64[],
        close=Float64[],
        high=Float64[],
        low=Float64[],
        timestamp=Date[],
        number_of_transactions=Int[]
    )

    # fill in the header dictionary -
    header_keys = [
        "ticker", "queryCount", "adjusted", "status", "request_id", "count"
    ]
    for key ∈ header_keys
        header_dictionary[key] = request_body_dictionary[key]
    end

    # populate the results DataFrame -
    results_array = request_body_dictionary["results"]
    for result_dictionary ∈ results_array
        
        # build a results tuple -
        result_tuple = (

            volume = result_dictionary["v"],
            volume_weighted_average_price = result_dictionary["vw"],
            open = result_dictionary["o"],
            close = result_dictionary["c"],
            high = result_dictionary["h"],
            low = result_dictionary["l"],
            timestamp = unix2datetime(result_dictionary["t"]*(1/1000)),
            number_of_transactions = result_dictionary["n"]
        )
    
        # push that tuple into the df -
        push!(df, result_tuple)
    end

    # return -
    return (header_dictionary, df)
end

function process_options_reference_call_response(body::String)

    # convert to JSON -
    request_body_dictionary = JSON.parse(body)

    # before we do anything - check: do we have an error?
    status_flag = request_body_dictionary["status"]
    if (status_flag == "ERROR")
        return polygon_error_handler(request_body_dictionary)
    end

    # initialize -
    header_dictionary = Dict{String,Any}()
    df = DataFrame(

        cfi=String[],
        contract_type=String[],
        exercise_style=String[],
        expiration_date=Date[],
        primary_exchange=String[],
        shares_per_contract=Int64[],
        strike_price=Float64[],
        ticker = String[],
        underlying_ticker = String[]
    )

    # fill in the header dictionary -
    header_keys = [
        "status", "request_id", "count", "next_url"
    ]
    for key ∈ header_keys
        header_dictionary[key] = request_body_dictionary[key]
    end

    # populate the results DataFrame -
    results_array = request_body_dictionary["results"]
    for result_dictionary ∈ results_array
        
        # build a results tuple -
        result_tuple = (

            cfi = result_dictionary["cfi"],
            contract_type = result_dictionary["contract_type"],
            exercise_style = result_dictionary["exercise_style"],
            expiration_date = Date(result_dictionary["expiration_date"]),
            primary_exchange = result_dictionary["primary_exchange"],
            shares_per_contract = result_dictionary["shares_per_contract"],
            strike_price = result_dictionary["strike_price"],
            ticker = result_dictionary["ticker"],
            underlying_ticker = result_dictionary["underlying_ticker"]
        )
    
        # push that tuple into the df -
        push!(df, result_tuple)
    end

    # return -
    return (header_dictionary, df)
end

function process_polygon_response(model::AbstractPolygonEndpointModel, 
    response::String)::Tuple

    # which handler should we call?
    if (isa(model,PolygonAggregatesEndpointModel) == true)
        return process_aggregates_polygon_call_response(response)
    elseif (isa(model, PolygonOptionsContractReferenceEndpoint) == true)
        return process_options_reference_call_response(response)
    end

    # default -
    return nothing
end

function polygon(base::String, model::AbstractPolygonEndpointModel; 
    handler::Function = process_polygon_response)::Tuple

    # build the url string -
    complete_url_string = build(base, model)

    # execute -
    result_string = http_get_call_with_url(complete_url_string) |> check

    # process and return -
    return handler(model, result_string)
end