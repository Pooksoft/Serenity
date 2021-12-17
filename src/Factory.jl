function _objective_function_sim(θ,X)

    # alias some parameters -
    α = θ[1]
    β = θ[2]
    R_firm = X[:,2]
    R_market = X[:,1]

    # compute the residual -
    ϵ = (R_firm .- α .- β*R_market)

    # return the total error -
    return sum(ϵ.^2)
end

function build_single_index_model(market::DataFrame, firm::DataFrame;
    risk_free_rate::Float64 = 5.15e-5)::SingleIndexModel

    # # market risk premimum -
    Y = (firm[!,:μ]) .- risk_free_rate

    # asset risk premimum -
    A = (market[!,:μ]) .- risk_free_rate

    # build the data array -
    number_of_time_steps = length(Y)
    X = [ones(number_of_time_steps,1) A]

    # compute the parameters -
    M = inv(transpose(X)*X)*transpose(X)
    theta_v = M*Y

	# compute the residuals -
    r = X*theta_v - Y
    μ_R = mean(r)
    σ_R = std(r)

	# build model -
	model = SingleIndexModel()
	model.α = theta_v[1]
	model.β = theta_v[2]
	model.r = risk_free_rate
	model.ϵ = Normal(μ_R, σ_R)

    # return -
    return model
end