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

function build_url_query_string(base::String, options::Dict{String,Any})::String

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

function process_aggregates_polygon_call_response(body::String)

    # convert to JSON -
    request_body_dictionary = JSON.parse(body)

    # initialize -
    header_dictionary = Dict{String,Any}()
    df = DataFrame(v=Float64[],vw=Float64[],o=Float64[],c=Float64[],h=Float64[],l=Float64[],t=Date[],n=Int[])
    #df = DataFrame()
    #df = Array{NamedTuple,1}()

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
            v = result_dictionary["v"],
            vw = result_dictionary["vw"],
            o = result_dictionary["o"],
            c = result_dictionary["c"],
            h = result_dictionary["h"],
            l = result_dictionary["l"],
            t = unix2datetime(result_dictionary["t"]*(1/1000)),
            n = result_dictionary["n"]
        )
    
        # push that tuple into the df -
        push!(df, result_tuple)
    end

    # return -
    return (header_dictionary, df)
end

function execute_polygon_aggregates_api_call(base::String, ticker::String, multiplier::Int,
    timespan::String, from::Date, to::Date, options::Dict{String,Any})

    # build up the base string -
    base_url = "$(base)/aggs/ticker/$(ticker)/range/$(multiplier)/$(timespan)/$(from)/$(to)?"

    # return -
    complete_url_string = build_url_query_string(base_url, options)

    # execute -
    result_string = http_get_call_with_url(complete_url_string) |> check

    # process and return -
    return (result_string |> process_aggregates_polygon_call_response)
end



