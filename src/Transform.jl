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