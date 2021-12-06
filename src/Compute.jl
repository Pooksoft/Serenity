function compute_return_array(data_table::DataFrame, key::Pair{Symbol,Symbol}; β = (1.0 / 365.0))


    # initialize -
    (number_of_rows, _) = size(data_table)
    return_table = DataFrame(timestep = Date[], P = Float64[], μ = Float64[])

    # main loop -
    for row_index = 2:number_of_rows

        # grab the date -
        tmp_date = data_table[row_index, key.first]

        # grab the price data -
        yesterday_close_price = data_table[row_index-1, key.second]
        today_close_price = data_table[row_index, key.second]

        # compute the diff -
        Δ = (1 / β) * log(today_close_price / yesterday_close_price)

        # push! -
        push!(return_table, (tmp_date, today_close_price, Δ))
    end

    # return -
    return return_table
end