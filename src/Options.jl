function description(asset::AbstractAsset)::String

    # get some data for this asset -
    ticker = asset.ticker
    sense = asset.sense

    if (isa(asset, CallContract) == true)
        
        strike_price = asset.K
        return "$(ticker) CALL@$(strike_price) $(sense)"
    
    elseif (isa(asset, PutContract) == true)

        strike_price = asset.K
        return "$(ticker) PUT@$(strike_price) $(sense)"
    
    elseif (isa(asset, Equity) == true)

        C = asset.C
        return "$(ticker) EQUITY@$(C) $(sense)"
    end
end

function payoff(contract::PutContract, 
    underlying_price::Float64)::Float64

    # get data from the contract -
    K = contract.K # contract strike price -

    # compute the intrinsic value for a put -
    return max(0.0, K - underlying_price)
end

function payoff(contract::CallContract, 
    underlying_price::Float64)::Float64

    # get data from the contract -
    K = contract.K # contract strike price -

    # compute the intrinsic value for a put -
    return max(0.0, underlying_price - K)
end

function payoff(contract::Equity, 
    underlying_price::Float64)::Float64

    # get data from this asset -
    C = contract.C

    # return -
    return (C - underlying_price)
end

function compute_profit_loss_at_expiration(assets::Array{AbstractAsset,1}, 
    underlying_prices::Array{Float64,1})::DataFrame

    # initialize -
    𝒜 = length(assets)
    𝒩 = length(underlying_prices)
    data_array = zeros(𝒩,(𝒜+2))
    col_name_array = Array{String,1}()
    
    # compute the profit and loss -
    for (asset_index, asset) ∈ enumerate(assets)
    
        # what is the sense of this asset?
        asset_sense = asset.sense
        C = asset.C  # cost of asset -


        for (price_index, price) ∈ enumerate(underlying_prices)

            # compute the intrinsic value of this contract -
            payoff_value = payoff(asset, price)

            plvalue = 0.0
            if asset_sense == :BUY
                plvalue = payoff_value - C # we purchased this asset/contract
            elseif asset_sense == :SELL
                plvalue = C - payoff_value # we sold this asset/contract
            end
            
            # capture -
            data_array[price_index,1] = price
            data_array[price_index, asset_index + 1] = plvalue
        end
    end

    for price_index ∈ 1:𝒩
        value = sum(data_array[price_index,2:end-1])
        data_array[price_index,end] = value
    end
    
    # construct the col names -
    push!(col_name_array,"P")
    for asset ∈ assets
        col_name = build_option_ticker_symbol(asset)
        push!(col_name_array, col_name)
    end
    push!(col_name_array,"Σ")

    # build the data frame -
    df = DataFrame(data_array,col_name_array)

    # return -
    return df
end