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

function compute_random_walk_model_trajectory(model::Distribution, initial_price::Float64,
    number_of_steps::Int64; number_of_sample_paths = 1)

    # initialize -
    number_of_steps = number_of_steps + 1
    price_array = Array{Float64,2}(undef, number_of_steps, number_of_sample_paths)

    # insert the first value -
    price_array[1, 1:number_of_sample_paths] .= initial_price

    # how many samples do we need?
    sample_return_array = rand(model, number_of_steps, number_of_sample_paths)

    # compute the price -
    for sample_path_index = 1:number_of_sample_paths
        for step_index = 2:number_of_steps
            price_array[step_index, sample_path_index] = price_array[step_index-1, sample_path_index] + sample_return_array[step_index, sample_path_index]
        end
    end

    # return -
    return price_array
end

function compute_rwm_cumulative_probabilty(compare::Function, price_array::Array{Float64,1})

    # initialize -
    number_of_samples = length(price_array)
    tmp_array = BitArray(undef, (number_of_samples, 1))

    # main -
    for sample_index = 1:number_of_samples

        # get the sample price -
        sample_price = price_array[sample_index]

        # check: which is larger, sample or target price?
        compare(sample_price) ? tmp_array[sample_index] = 1 : tmp_array[sample_index] = 0
    end

    # sum the tmp_array -
    number_of_larger_values = sum(tmp_array)

    # compute the probability -
    return (number_of_larger_values / number_of_samples)
end

function compute_rwm_cumulative_probabilty(price_array::Array{Float64,1}, target_price::Float64)

    # initialize -
    number_of_samples = length(price_array)
    tmp_array = Array{Int64,1}()

    # main -
    for sample_index = 1:number_of_samples

        # get the sample price -
        sample_price = price_array[sample_index]

        # check: which is larger, sample or target price?
        sample_price <= target_price ? push!(tmp_array, 1) : push!(tmp_array, 0)
    end

    # sum the tmp_array -
    number_of_larger_values = sum(tmp_array)

    # compute the probability -
    return (number_of_larger_values / number_of_samples)
end

# short cut method RWMC method -
(model::Distribution)(initial_price::Float64, number_of_steps::Int64; number_of_sample_paths = 1) =
    compute_random_walk_model_trajectory(model, initial_price, number_of_steps;
        number_of_sample_paths = number_of_sample_paths)

