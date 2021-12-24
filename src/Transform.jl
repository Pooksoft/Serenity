function transform_date_fields(data_frame::DataFrame, col_key::Symbol)

    # creates a col of string dates -
    date_col = string.(data_frame[!, col_key])
    tmp_array = Array{Date,1}()
    for (_, date_value) in enumerate(date_col)

        # split around "/"
        date_components = split(date_value, "/")

        # process the month -
        month_value = date_components[1]
        if (length(month_value) < 2)
            month_value = lpad(month_value, 1, "0")
        end

        # process the day -
        day_value = date_components[2]
        if (length(day_value) < 2)
            day_value = lpad(day_value, 1, "0")
        end

        # process the year -
        year_value = date_components[3]
        if (length(year_value) < 3)
            year_value = "20$(year_value)"
        end

        # new date -
        new_date_string = "$(year_value)-$(month_value)-$(day_value)"
        new_date_object = Dates.Date(new_date_string)
        push!(tmp_array, new_date_object)
    end

    # remove the col_key -
    select!(data_frame, Not(col_key))
    insertcols!(data_frame, col_key => tmp_array)
end

function filter(data::DataFrame, start_date::Date, end_date::Date;
    key::Symbol = :timestamp)::DataFrame

    return filter(key => (x) -> (x>=start_date && x<=end_date), data);
end

function extract_data_block_for_date_range(data_table::DataFrame, start_date::Date, end_date::Date;
    key::Symbol = :timestamp)::DataFrame

    # find the range -
    idx_range = findall(x->(x>=start_date && x<=end_date), data_table[!, key]) |> collect

    # extract the range -
    data_block = data_table[idx_range,:]

    # return -
    return data_block
end

function reverse_row_order_in_table(table::DataFrame)

    # how many rows do we have?
    (number_of_rows, _) = size(table);

    # create a data array -
    tmp_date_array = Array{Date,1}()
    for row_index ∈ 1:number_of_rows
    end

    df_tmp = DataFrame()
    for row_index ∈ number_of_rows:-1:1
        df_row = table[row_index,:]
        push!(df_tmp, df_row)
    end

    return df_tmp
end