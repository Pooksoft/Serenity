mutable struct SingleIndexModel

    # model -
    α::Float64          # firm specific unexplained return
    β::Float64          # relationship between the firm and the market
    r::Float64          # risk free rate of return 
    ϵ::Distribution     # random shocks 

    # constructor -
    SingleIndexModel() = new()
end

mutable struct GeometricBrownianMotionModel

    # model parameters -
    μ::Float64
    σ::Float64

    GeometricBrownianMotionModel() = new()
end

# --- ASSETS ------------------------------------------------- #
abstract type AbstractAsset end
mutable struct CallContract <: AbstractAsset

    # data -
    ticker::String
    expiration::Date
    K::Float64
    C::Float64
    sense::Symbol
    
    # constructor
    CallContract() = new()
end

mutable struct PutContract <: AbstractAsset

    # data -
    ticker::String
    expiration::Date
    K::Float64
    C::Float64
    sense::Symbol
    
    # constructor
    PutContract() = new()
end

mutable struct Equity <: AbstractAsset

    # data -
    ticker::String
    C::Float64
    sense::Symbol
    
    # constructor
    Equity() = new()
end



# --- Polygon API model --------------------------------------------- #
abstract type AbstractPolygonEndpointModel end
mutable struct PolygonAggregatesEndpointModel <: AbstractPolygonEndpointModel

    # data -
    adjusted::Bool
    limit::Int64
    sortdirection::String
    apikey::String
    ticker::String
    to::Date
    from::Date
    multiplier::Int64
    timespan::String

    # constructor -
    PolygonAggregatesEndpointModel() = new()
end

# Setup call for options data -
@enum PolygonOptionSortKey ticker=1 underlying_ticker expiration_date strike_price
@enum PolygonOptionContractType call=1 put 
mutable struct PolygonOptionsContractReferenceEndpoint <: AbstractPolygonEndpointModel

    # data -
    ticker::Union{Nothing,String}
    underlying_ticker::Union{Nothing, String}
    contract_type::PolygonOptionContractType
    expiration_date::Date
    limit::Int64
    order::String
    sort::PolygonOptionSortKey
    apikey::String

    # constructor -
    PolygonOptionsContractReferenceEndpoint() = new()
end