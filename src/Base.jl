function build_url_query_string(base::String, options::Dict{String,String})::String

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

function check(result::Some)::(Union{Nothing,T} where {T<:Any})

    # ok, so check, do we have an error object?
    # Yes: log the error if we have a logger, then throw the error. 
    # No: return the result.value

    # Error case -
    if (isa(something(result), Exception) == true)

        # get the error object -
        error_object = result.value

        # get the error message as a String -
        error_message = sprint(showerror, error_object, catch_backtrace())
        @error(error_message)

        # throw -
        throw(result.value)
    end

    # default -
    return result.value
end

function process_csv_api_data(api_call_raw_data::String)::DataFrame

    # create a data table from the CSV data -
    tmp_data_table = CSV.read(IOBuffer(api_call_raw_data), DataFrame)

    # sort the table according to the timestamps -
    idx_sort = sortperm(tmp_data_table[:, 1])

    # create a sorted data table -
    sorted_data_table = tmp_data_table[idx_sort, :]

    # return -
    return sorted_data_table
end