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