### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f9e5092b-a158-48cb-a919-96f30ec51038
begin

    # setup: we need to specify the paths where we can find our project resources
    _PATH_TO_ROOT = pwd()
    _PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")
    _PATH_TO_DATA = joinpath(_PATH_TO_ROOT, "data")
    _PATH_TO_CONFIG = joinpath(_PATH_TO_ROOT, "configuration")

    # what packages are we going to use?
    using DataFrames
    using CSV
    using HTTP
    using TOML
    using ProgressMeter
    using Dates
    using StatsPlots
    using Statistics
    using Distributions
    using PlutoUI
    using PlotThemes
	using HypothesisTests
	using PrettyTables
	
    theme(:default)

    # load my lib of code -
    include(joinpath(_PATH_TO_SRC, "Include.jl"))

    # alias the AlphaVantage URL -
    DATASTORE_URL_STRING = "https://www.alphavantage.co/query?"

    # What we have here is a classic good news, bad news situation ...
    # Bad news: We don't check in our AlphaVantage API key to GitHub (sorry). 
    # Good news: AlphaVantage API keys are free. 
    # check out: https://www.alphavantage.co/support/#api-key
    configuration_dictionary = TOML.parsefile(joinpath(_PATH_TO_CONFIG, "Configuration.toml"))

    # get some stuff from the configuration dictionary -
    ALPHAVANTAGE_API_KEY = configuration_dictionary["API"]["alphavantage_api_key"]

    # background color plots -
    background_color_outside = RGB(1.0, 1.0, 1.0)
    background_color = RGB(0.99, 0.98, 0.96)
    CB_BLUE = RGB(68 / 255, 119 / 255, 170 / 255)
    CB_LBLUE = RGB(102 / 255, 204 / 255, 238 / 255)
    CB_GRAY = RGB(187 / 255, 187 / 255, 187 / 255)
    CB_RED = RGB(238 / 255, 102 / 255, 119 / 255)

    # show -
    nothing
end

# ╔═╡ 2bb52ee4-1c6f-46b6-b105-86827ada0f75
md"""
# Random Walk Simulation of Risky Asset Prices 
"""

# ╔═╡ 2f499c95-38cf-4856-b199-6c9aac44237a
md"""
### Introduction
"""

# ╔═╡ 34bf07e8-5c47-4aa8-aa2c-161709c158be
md"""
### Materials and Methods

##### Random Walk Models (RWMs)
[Random walk models](https://en.wikipedia.org/wiki/Random_walk) are tools for computing the time evolution of the prices of risky assets, for example, the price of a stock with the ticker symbol `XYZ`. Let the price of `XYZ` at time $j$ be given by $P_{j}$ (units: USD/share). Then during the next time period (index $j+1$) the `XYZ` share price will be given by:

$$P_{j+1} = P_{j}\exp\left(\mu_{j\rightarrow{j+1}}\Delta{t}\right)$$

where $\Delta{t}$ denotes the length of time between periods $j$ and $j+1$ (units: time) and $\mu_{j\rightarrow{j+1}}$ denotes the growth rate (in the financial world called the [return](https://www.investopedia.com/terms/a/absolutereturn.asp)) between time period(s) $j$ and $j+1$ (units: time$^{-1}$). The growth rate $\mu_{j\rightarrow{j+1}}$ is different for different tickers, and is not constant between time period $j\rightarrow{j+1}$. Dividing both sides by $P_{j}$ and taking the [natural log](https://en.wikipedia.org/wiki/Natural_logarithm) gives the expression:

$$p_{j+1} = p_{j} + \mu_{j\rightarrow{j+1}}\Delta{t}$$

where $p_{\star}$ denotes the log of the price at time $\star$ (units: dimensionless). This expression is a basic [random walk model](https://en.wikipedia.org/wiki/Random_walk). The challenge to using this model is how to estimate the growth rate 
$\mu_{j\rightarrow{j+1}}$ because the values for $\mu_{j\rightarrow{j+1}}$ are [time-dependent random variables](https://en.wikipedia.org/wiki/Stochastic_process) which are [distributed](https://en.wikipedia.org/wiki/Probability_distribution) in some unknown way. However, if we had a value for $\mu_{j\rightarrow{j+1}}$ (or at least a model to estimate a value for it, denoted as $\hat{\mu}_{j\rightarrow{j+1}}$), we could 
simulate the price of `XYZ` during the next time period $j+1$ given knowledge of the price now (time period $j$). For example, given a value for the close price of `XYZ` on Monday, we could estimate a value for the close price of `XYZ` on Tuesday if we had a model for the daily growth rate.

##### Estimating models for the growth rate $\mu$ from historical data
Suppose values the growth rate were governed by some [probability distribution](https://en.wikipedia.org/wiki/Probability_distribution) $\mathcal{D}\left(\bar{m},\sigma\right)$ where $\bar{m}$ 
denotes the [mean value](https://en.wikipedia.org/wiki/Mean) and $\sigma$ denotes the [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) of the growth rate $\mu_{j\rightarrow{j+1}}$. Of course, we do not know the actual values for 
$\left(\bar{m},\sigma\right)$, nor do we know the form of $\mathcal{D}\left(\bar{m},\sigma\right)$, but we can estimate them from historial price data. An even better way to learn $\mathcal{D}\left(\bar{m},\sigma\right)$ would be to estimate the form of the distribution from the data itself using a technique such as [Kernel Density Estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation) or KDE. However, KDE is a little advanced for now, maybe later. 

"""

# ╔═╡ 9943d000-83d0-413d-a231-0295fb19df71
md"""
### Results and Discussion
"""

# ╔═╡ f66a480b-3f0c-4ebf-a8b8-e0f91dff851d
md"""
##### Download historical price data from the [Alphavantage.co](https://www.alphavantage.co) financial data application programming interface (API)
We specify a list of ticker symbols that we want to model. Next, we check to see if we have price data already saved locally for each ticker. 
* If yes, then we load the saved file as a [DataFrame](https://dataframes.juliadata.org/stable/) and store it in the `price_data_dictionary` where the keys are the ticker strings e.g., MSFT, etc. 
* if no, we download new data from [Alphavantage.co](https://www.alphavantage.co), save this data locally as a `CSV` file (<ticker>.csv), and finally we store the price data as a DataFrame in the `price_data_dictionary` where the keys are the ticker symbols.
The parameters for the [Alphavantage.co](https://www.alphavantage.co) API call are stored in the `query_parameters` dictionary. For estimating a SIM, we use the daily close price for the last 100 trading days for each ticker. 
"""

# ╔═╡ 54efa70c-bac6-4d7c-93df-0dfd1b89769d
# Pooksoft Industrial Average (PSIA) -> the DJIA + some stuff
ticker_symbol_array = sort(["MSFT", "ALLY", "MET", "AAPL", "GM", "PFE", "TGT", "WFC", "AIG", "F", "GE", "AMD",
    "MMM", "AXP", "AMGN", "BA", "CAT", "CVX", "CSCO", "KO", "DIS", "DOW", "GS", "HD", "IBM", "HON", "INTC", "JNJ", "JPM",
    "MCD", "MRK", "NKE", "PG", "CRM", "TRV", "UNH", "VZ", "V", "WBA", "WMT"
]);

# ╔═╡ 61ab2949-d72f-4d80-a717-4b6a9227de0e
md"""
##### Estimate the growth rate distributions from historical data

To estimate a model for the return $\mathcal{D}\left(\bar{m},\sigma\right)$, we first compute the historical return for a ticker. 
The [Serenity library](https://github.com/Pooksoft/Serenity.git) provides a method to compute the [log return](https://en.wikipedia.org/wiki/Rate_of_return) from historical data (contained in a DataFrame):

	compute_log_return_array(data::DataFrame,map::Pair{Symbol,Symbol}; Δt::Float64 = 1.0) -> DataFrame

	Computes μ of historical price data. 
	
	Arguments:
	data  	DataFrame that holds historical price data where each row is a record
	map 	Different APIs return data with different field names. The map arg connects the time and price fields
	Δt 		Time step size. Default is 1.0 in the units you want the μ to be calculated in
"""

# ╔═╡ a39b90ec-a4c5-472f-b00e-c81ee9c5576f
md"""
__Fig. 1__: Histogram of the adjusted daily returns computed using the `stephist` routine of the [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl) package for the 𝒫 = 40 ticker symbols in the PSIA. The number of bins was set to be 10% of the number of historical records for `AAPL` (shown in red) and was constant for each ticker. The returns were calculated using approximately 10 to 20 years of historical daily adjusted close price data (depending upon the ticker symbol) 
downloaded from using [Alphavantage.co](https://www.alphavantage.co) application programming interface. 
"""

# ╔═╡ edfbf364-e126-4e95-93d2-a6adfb340045
md"""
###### The daily return of individual ticker symbols is not normally distributed

A key assumption often made in mathematical finance is that the growth rate $\mu_{j\rightarrow{j+1}}$ (return) computed using the log of the prices is normally distributed. However, a closer examination of the $\mu_{j\rightarrow{j+1}}$ values computed from historical data suggests this is not true (Fig. 2). Visual inspection of the histograms computed using a [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution) (Fig 2, dark blue) for the returns appear to be closer to the historical data (Fig 2, red) compared to a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) (Fig 2, light blue). 
"""

# ╔═╡ a786ca10-06d2-4b76-97a9-2bcf879ea6cb
# fit a distribution to a ticker -
single_asset_ticker_symbol = "GS";

# ╔═╡ a3d29aa3-96ca-4681-960c-3b4b04b1e40d
md"""
__Fig 2__: Comparison of the actual (red line) and estimated histogram for the adjusted daily returns for ticker = $(single_asset_ticker_symbol) using a [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution) for the return model (blue line) and a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) for the return model (light blue line). The return model distributions were estimated using the `fit` routine of the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package. The `fit` routine uses multiple approaches to estimate the parameters in the distribution including [Maximum Likelihood Estimation (MLE)](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation). 
"""

# ╔═╡ 1d72b291-24b7-4ec6-8307-1da0bc4a9183
md"""
To further explore the question of a Normal versus Laplace return distribution, we constructed [QQ-Plots](https://en.wikipedia.org/wiki/Q–Q_plot) to visualize the historical return data versus a Laplace distribution (Fig. 3), and performed the [Kolmogorov–Smirnov (KS) test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test), using the KS test implementation in the [HypothesisTests.jl](https://juliastats.org/HypothesisTests.jl/latest/nonparametric/#Kolmogorov-Smirnov-test-1) package, to test whether the historical return data were Normally or Laplace distributed (Table 1). 
"""

# ╔═╡ 1867ecec-3c3d-4b2b-9036-1488e2184c40
md"""
__Fig 3__: Comparison of the actual (red line) and estimated histogram for the adjusted daily returns for ticker = $(single_asset_ticker_symbol) using a [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution) for the return model (blue line) and a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) for the return model (light blue line). The return model distributions were estimated using the `fit` routine of the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package. The `fit` routine uses multiple approaches to estimate the parameters in the distribution including [Maximum Likelihood Estimation (MLE)](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation). 
"""

# ╔═╡ 379373b1-e563-4341-9978-5b35c768b5c7
md"""
__Table 1__: [Kolmogorov–Smirnov (KS)](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test) test results for Normal and Laplace distribution for stocks in the PSIA (𝒫 = 40). The historical return values for all tickers in PSIA were not Normally distributed. However, the returns for 23 of the 40 members of the PSIA were governed by a Laplace distribution. On the other hand, 17 tickers failed the KS test. Thus, the returns for these tickers were not governed by a Normal or Laplace distributions.   
"""

# ╔═╡ 10fa507e-1429-4eb0-b74c-1e6638725690
md"""
##### Monte Carlo simulations of the daily close price using a random walk model

A basic explain of Monte Carlo simulations goes here
"""

# ╔═╡ e36979d5-c1b6-4c17-a65a-d8de8e6bd8d0
md"""
Show actual price trajectory? $(@bind show_real_traj CheckBox()) 
"""

# ╔═╡ c32725a4-e276-4372-8d06-d40ba52c9f09
md"""
### Conclusions


"""

# ╔═╡ f1a71f47-fb19-4988-a439-2ff8d38be5b7


# ╔═╡ a6c4e663-f1e3-4e0c-a8bf-7c13fcb076f0
function download_ticker_data(ticker_symbol_array::Array{String,1})::Dict{String,DataFrame}

    # initialize some storage -
    price_data_dictionary = Dict{String,DataFrame}() # Dict key => value where key = ticker symbol and value = price DataFrame

    # setup: specify the query parameters for the call to Alphavantage -
    # check out the API documentation: https://www.alphavantage.co/documentation/
    query_parameters = Dict{String,String}()
    query_parameters["function"] = "TIME_SERIES_DAILY_ADJUSTED" # get the daily close price adjusted
    query_parameters["apikey"] = ALPHAVANTAGE_API_KEY # your alphavantage API key goes here (see setup/config block below)
    query_parameters["datatype"] = "csv" # download in CSV (the other option is JSON) format
    query_parameters["outputsize"] = "full" # compact -vs- full: last 100 trading days -vs- up to 20 years of data

    # data type: sub dir where we save the data
    data_type_flag = "daily"
    training_prediction_flag = "training"

    # main download loop -
    for ticker_symbol in ticker_symbol_array

        # check -> have we downloaded this data already?
        data_file_path = joinpath(_PATH_TO_DATA, "$(training_prediction_flag)", "$(data_type_flag)",
            "$(ticker_symbol).csv")
        if (ispath(data_file_path) == false || does_data_file_exist(data_file_path) == false)

            # let use know what is going on ...
            with_terminal() do
                println("Starting $(ticker_symbol) download ...")
            end

            # We do NOT have this data file -> download from Alphavantage.co
            # execute API call -> check to see if error -> turn into DataFrame
            query_parameters["symbol"] = ticker_symbol
            url_string = build_url_query_string(DATASTORE_URL_STRING, query_parameters)
            data_table = http_get_call_with_url(url_string) |> check |> process_csv_api_data

            # sleep -
            sleep(10) # to avoid issue w/API limits

            # make the dir to save the data (if we don't have it already) -
            mkpath(joinpath(_PATH_TO_DATA, "$(training_prediction_flag)", "$(data_type_flag)"))

            # dump data table to disk -
            CSV.write(data_file_path, data_table)

            # put the price DataFrame into the price_data_dictionary -
            price_data_dictionary[ticker_symbol] = data_table

        else

            # we already have this ticker downloaded. load the existing file from disk -
            data_table = CSV.read(data_file_path, DataFrame)

            # put the price DataFrame into the price_data_dictionary -
            price_data_dictionary[ticker_symbol] = data_table
        end
    end

    # show -
    return price_data_dictionary
end

# ╔═╡ 5c5d5eeb-6775-452f-880d-7b4fa2acda57
# download the data (or load from local cache)
price_data_dictionary = download_ticker_data(ticker_symbol_array);

# ╔═╡ 34b06415-21c1-4904-97f0-ab614447355c
begin

    # initialize -
    return_data_dictionary = Dict{String,DataFrame}()

    # compute -
    for ticker_symbol ∈ ticker_symbol_array

        # compute_return_array function is provided by Serenity -> computes the log return given a DataFrame, returns a DataFrame
        return_data_dictionary[ticker_symbol] = compute_log_return_array(price_data_dictionary[ticker_symbol],
            :timestamp => :adjusted_close; Δt = (1.0 / 1.0))
    end
end

# ╔═╡ cbbd8670-49ab-4601-b8d7-9f3f456752e8
begin

    # compute some ranges -
    ticker_symbol = "AAPL"
    μ_avg = mean(return_data_dictionary[ticker_symbol][!, :μ])

    # number of bins - 10% of the length of a test ticker
    number_of_bins = convert(Int64, (floor(0.1 * length(return_data_dictionary[ticker_symbol][!, :μ]))))

    # number of tickers -
    number_of_ticker_symbols = length(ticker_symbol_array)
    for ticker_index = number_of_ticker_symbols:-1:1
        local_ticker_symbol = ticker_symbol_array[ticker_index]

        if (ticker_index == number_of_ticker_symbols)

            # make the first plot -
            stephist(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true,
                background_color = background_color, background_color_outside = background_color_outside,
                foreground_color_minor_grid = RGB(1.0, 1.0, 1.0),
                lw = 1, c = CB_GRAY, foreground_color_legend = nothing, label = "")

        elseif (ticker_index != 1 && ticker_index != number_of_ticker_symbols)

            stephist!(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true, lw = 1,
                c = CB_GRAY, label = "")
        else
            stephist!(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true, lw = 2,
                c = CB_RED, label = "$(local_ticker_symbol)")
        end
    end

    # label the plots -
    xlabel!("Daily return μ (1/day)", fontsize = 18)
    ylabel!("Frequency (N=$(number_of_bins); dimensionless)", fontsize = 14)
    xlims!((-100.0 * μ_avg, 100.0 * μ_avg))
end

# ╔═╡ 6bf06c12-cf25-43c4-81f3-b1d79d13fc94
begin

    # get the historical return data -
    μ_vector = return_data_dictionary[single_asset_ticker_symbol][!, :μ]

    # fit the model -
    NORMAL = fit(Normal, μ_vector)
    LAPLACE = fit(Laplace, μ_vector)

    # generate 10_000 samples
    S = rand(LAPLACE, 20000)
    SN = rand(NORMAL, 20000)

    # plot against actual -
    stephist(return_data_dictionary[single_asset_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true, 
		lw = 2, c = CB_RED,
        label = "$(single_asset_ticker_symbol)", background_color = background_color,
        background_color_outside = background_color_outside, foreground_color_legend = nothing)

    stephist!(S, bins = number_of_bins, normed = :true, lw = 2, c = CB_BLUE,
        label = "$(single_asset_ticker_symbol) Laplace Model")

    stephist!(SN, bins = number_of_bins, normed = :true, lw = 2, c = CB_LBLUE,
        label = "$(single_asset_ticker_symbol) Normal Model")

    xlabel!("Daily return μ (1/day)", fontsize = 18)
    ylabel!("Frequency (N=$(number_of_bins); dimensionless)", fontsize = 14)
    xlims!((-100.0 * μ_avg, 100.0 * μ_avg))
end

# ╔═╡ 0e09d312-2ddf-4d1f-8ad5-a50fb48ca4dd
begin
	qqplot(Laplace, μ_vector, label="QQPlot $(single_asset_ticker_symbol)", legend=:topleft, 
		foreground_color_legend = nothing, c=CB_LBLUE, markerstrokecolor=CB_BLUE, lw=2, 
		background_color = background_color, background_color_outside = background_color_outside)
	xlabel!("Theoretical Quantiles", fontsize=18)
	ylabel!("Sample Quantiles", fontsize=18)
end

# ╔═╡ f0ee3633-1d28-4344-9b85-f5433679582a
test_result = HypothesisTests.ExactOneSampleKSTest(unique(μ_vector), LAPLACE)

# ╔═╡ 8b2725a2-8007-46b8-a160-75042562794d
with_terminal() do

	number_of_ticker_symbols = length(ticker_symbol_array)
	state_table = Array{Any,2}(undef, number_of_ticker_symbols, 5)
	for ticker_index ∈ 1:number_of_ticker_symbols

		my_ticker_symbol = ticker_symbol_array[ticker_index]
		μ_vector = return_data_dictionary[my_ticker_symbol][!, :μ]
		LAPLACE = fit(Laplace, μ_vector)
		NORMAL = fit(Normal, μ_vector)
		test_result_laplace = HypothesisTests.ExactOneSampleKSTest(unique(μ_vector), LAPLACE)
		test_result_normal = HypothesisTests.ExactOneSampleKSTest(unique(μ_vector), NORMAL)
		
		state_table[ticker_index,1] = my_ticker_symbol
		state_table[ticker_index,2] = pvalue(test_result_normal)
		state_table[ticker_index,3] = pvalue(test_result_normal)≥0.05 ? true : false
		state_table[ticker_index,4] = pvalue(test_result_laplace)
		state_table[ticker_index,5] = pvalue(test_result_laplace)≥0.05 ? true : false
	end

	table_header = (
		["ticker", "Normal", "pₙ ≥ 0.05", "Laplace", "pₗ ≥ 0.05"],
		["", "p-value", "95% CI", "p-value", "95% CI"]
	)
	pretty_table(state_table, header=table_header)
	
end

# ╔═╡ 5a3500c2-4f82-43e9-a31b-d530f56fdbe9
begin
    
	# how many steps, sample paths etc -
    number_of_days = 21
    number_of_sample_paths = 25000
	monte_carlo_simulation_dictionary = Dict{String,Array{Float64,2}}()

	for ticker_symbol ∈ ticker_symbol_array

		# what is the *actual* price data?
    	actual_price_data = price_data_dictionary[ticker_symbol][end-number_of_days:end, :adjusted_close]

    	# get initial price -
    	initial_price_value = log(actual_price_data[1])

    	# compute a set of possible trajectories -> convert back to actual price -
    	monte_carlo_simulation_dictionary[ticker_symbol] = LAPLACE(initial_price_value, number_of_days;
        		number_of_sample_paths = number_of_sample_paths) .|> exp
	end
    
    # show -
    nothing
end

# ╔═╡ a1e1d5f8-e06e-4682-ab54-a9454a8e3b30
md"""
__Fig 4__: In sample random walk simulation of ticker = $(single_asset_ticker_symbol) for 𝒯 = $(number_of_days) days. Blue lines denotes simulated sample paths while the red line denotes actual price trajectory for ticker $(single_asset_ticker_symbol). The simulation consisted on N = $(number_of_sample_paths) sample paths.
"""

# ╔═╡ b547311c-ddf0-4053-9de4-f0e85b861e63
md"""
__Table 2__: Comparison of the actual close price versus the Monte Carlo simulated close price for a 𝒯 = $(number_of_days) day prediction horizon for each ticker in the PSIA (𝒫 = 40). Each ticker was classified into class c ∈ {-1,0,1} where: +1 overbought, 0 in-range, or -1 oversold. The classification was based upon whether the actual close price Pₐ ∈ Pₑ ± σ, where Pₐ denotes the actual close price (units: USD/share), Pₑ denotes the mean simulated close price (units: USD/share), and σ denotes the standard deviation of the simulated close price (units: USD/share) computed over the family Monte Carlo trajectories (N = $(number_of_sample_paths)).
"""

# ╔═╡ 849f69b0-07af-40ab-8295-c0b80a26a2d5
md"""
__Fig 5__: In sample random walk simulation (N = $(number_of_sample_paths) sample paths) of ticker $(single_asset_ticker_symbol) for T = $(number_of_days) days. Gray lines denotes simulated sample paths. Red line denotes actual price trajectory for ticker $(single_asset_ticker_symbol).
"""

# ╔═╡ aeafe1ed-f217-48fd-9624-add5f6f791e6
begin
	
	# setup some preliminaries -
	skip_factor = convert(Int64, (floor(0.005 * number_of_sample_paths)))
    plot_index_array = 1:skip_factor:number_of_sample_paths |> collect
	simulated_price_trajectory = monte_carlo_simulation_dictionary[single_asset_ticker_symbol]
	
	# what is the *actual* price data?
    actual_price_data = price_data_dictionary[single_asset_ticker_symbol][end-number_of_days:end, :adjusted_close]
	
    # plot -
    plot(simulated_price_trajectory[:, plot_index_array], c = CB_LBLUE, legend = false, label = "", lw = 1,
        background_color = background_color, background_color_outside = background_color_outside)

    if (show_real_traj == true)
        plot!(actual_price_data[1:end], c = CB_RED, lw = 3, legend = :topleft, 
			label = "$(single_asset_ticker_symbol) actual", foreground_color_legend = nothing)
    end

    xlabel!("Time step index (day)", fontsize = 18)
    ylabel!("Simulated $(single_asset_ticker_symbol) close price (USD/share)", fontsize = 14)
    # title!("Random walk simulation $(single_asset_ticker_symbol) (N = $(number_of_sample_paths))", fontsize=12)
end

# ╔═╡ 36372d31-215d-4299-b4f1-49e42d8b0dbd
begin

	estimated_mean_price = round(mean(simulated_price_trajectory[end, :]), sigdigits = 4)
    std_estimated_price = round(std(simulated_price_trajectory[end, :]), sigdigits = 4)
	
    # compute the cumulative probability in the range [0,μ + 3*σ]
    LB = estimated_mean_price - 3 * std_estimated_price
    UB = estimated_mean_price + 3 * std_estimated_price
    price_range = range(LB, stop = UB, length = 1000) |> collect
    cprob = Array{Float64,1}()

    for price in price_range
        p = compute_rwm_cumulative_probabilty(x -> (x <= price), simulated_price_trajectory[end, :])
        push!(cprob, p)
    end

    # plot -
    #plot(price_range, cprob, legend=:right, label="P(X≤x)", lw=2)
    plot(price_range, 1 .- cprob, legend = :topright, label = "P(X≥x)", lw = 2, c = :red,
        background_color = background_color, background_color_outside = background_color_outside, 
		foreground_color_legend = nothing)
    xlabel!("$(single_asset_ticker_symbol) close daily price (USD/share)", fontsize = 18)
    ylabel!("1 - cumulative probability P(X≤x)")
    #title!("$(single_asset_ticker_symbol) (T = $(number_of_days) days)", fontsize=18)
end

# ╔═╡ 2ed2f3c3-619b-4aed-b88b-b92a43578d84
with_terminal() do

	# how many tickers do we have?
	number_of_ticker_symbols = length(ticker_symbol_array)

	# initialize some storage -
	state_table = Array{Any,2}(undef, number_of_ticker_symbols, 8)
	
	for ticker_symbol_index ∈ 1:number_of_ticker_symbols 

		# get the symbol -
		ticker_symbol = ticker_symbol_array[ticker_symbol_index]

		# compute some data re this symbol -
		simulated_price_trajectory = monte_carlo_simulation_dictionary[ticker_symbol]
		estimated_mean_price = round(mean(simulated_price_trajectory[end, :]), sigdigits = 4)
    	std_estimated_price = round(std(simulated_price_trajectory[end, :]), sigdigits = 4)
		actual_price_data = price_data_dictionary[ticker_symbol][end-number_of_days:end, :adjusted_close]
    	price_actual = actual_price_data[end]
		LB = estimated_mean_price - std_estimated_price
        UB = estimated_mean_price + std_estimated_price
		
		# populate state table -
		state_table[ticker_symbol_index,1] = ticker_symbol
		state_table[ticker_symbol_index,2] = price_actual
		state_table[ticker_symbol_index,3] = estimated_mean_price
		state_table[ticker_symbol_index,4] = LB
		state_table[ticker_symbol_index,5] = UB
		
        δL = max(0,round(LB - price_actual, sigdigits = 4))
        δU = max(0,round(price_actual - UB, sigdigits = 4))
		state_table[ticker_symbol_index,6] = δL
		state_table[ticker_symbol_index,7] = δU	

		if (δL == 0.0 && δU == 0.0)
			state_table[ticker_symbol_index,8] = 0 	# exepcted
		elseif (δL>0.0 && δU == 0.0)
			state_table[ticker_symbol_index,8] = -1 # oversold
		elseif (δL == 0.0 && δU > 0.0)
			state_table[ticker_symbol_index,8] = 1 # overbought
		end
	end

	table_header = (
		["ticker", "actual Pₐ", "estimate Pₑ", "ℒ = (Pₑ - σ)", "𝒰 = (Pₑ + σ)", "δL = max(0, ℒ - Pₐ)", 
			"δU = max(0, Pₐ - 𝒰)", "class c"],
		["", "USD/share", "USD/share", "USD/share", "", "USD/share", "USD/share","c ∈ {-1,0,1}"]
	)
	
	pretty_table(state_table, header=table_header)
	
end

# ╔═╡ c8c3fe32-560d-11ec-0617-2dc33608384a
html"""
<style>
main {
    max-width: 1200px;
    width: 85%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}

.H1 {
    padding: 0px 30px;
}
</style>"""

# ╔═╡ 5394b13a-e629-47c8-902d-685c061b37ae
html"""
<script>
	// initialize -
	var section = 0;
	var subsection = 0;
	var headers = document.querySelectorAll('h3, h5');
	
	// main loop -
	for (var i=0; i < headers.length; i++) {
	    
		var header = headers[i];
	    var text = header.innerText;
	    var original = header.getAttribute("text-original");
	    if (original === null) {
	        
			// Save original header text
	        header.setAttribute("text-original", text);
	    } else {
	        
			// Replace with original text before adding section number
	        text = header.getAttribute("text-original");
	    }
	
	    var numbering = "";
	    switch (header.tagName) {
	        case 'H3':
	            section += 1;
	            numbering = section + ".";
	            subsection = 0;
	            break;
	        case 'H5':
	            subsection += 1;
	            numbering = section + "." + subsection;
	            break;
	    }
		// update the header text 
		header.innerText = numbering + " " + text;
	};
</script>"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[compat]
CSV = "~0.9.11"
DataFrames = "~1.3.0"
Distributions = "~0.25.35"
HTTP = "~0.9.17"
HypothesisTests = "~0.10.6"
PlotThemes = "~2.0.1"
PlutoUI = "~0.7.22"
PrettyTables = "~1.2.3"
ProgressMeter = "~1.7.1"
StatsPlots = "~0.14.29"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra"]
git-tree-sha1 = "2ff92b71ba1747c5fdd541f8fc87736d82f40ec9"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.4.0"

[[deps.Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "49f14b6c56a2da47608fe30aed711b5882264d7a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.11"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "2e993336a3f68216be91eb8ee4625ebbaba19147"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "d6cc7abd52ebae5815fd75f6004a44abcf7a6b00"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.35"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "74ef6288d071f58033d54fd6708d4bc23a8b8972"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Random", "Rmath", "Roots", "Statistics", "StatsBase"]
git-tree-sha1 = "dc9bb7abfa265e0cf030635315184a476a2dd5f3"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.10.6"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "ca99cac337f8e0561c6a6edeeae5bf6966a78d21"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "8d958ff1854b166003238fe191ec34b9d592860a"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.8.0"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "65ebc27d8c00c84276f14aaf4ff63cbe12016c70"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "565564f615ba8c4e4f40f5d29784aa50a8f7bbaf"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.22"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "51ee572776905ee34c0568f5efe035d44bf59f74"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.11"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "bedb3e17cc1d94ce0e6e66d3afa47157978ba404"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.14"

[[deps.StatsPlots]]
deps = ["Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "d6956cefe3766a8eb5caae9226118bb0ac61c8ac"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.14.29"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "80661f59d28714632132c73779f8becc19a113f2"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─2bb52ee4-1c6f-46b6-b105-86827ada0f75
# ╟─2f499c95-38cf-4856-b199-6c9aac44237a
# ╟─34bf07e8-5c47-4aa8-aa2c-161709c158be
# ╟─9943d000-83d0-413d-a231-0295fb19df71
# ╟─f66a480b-3f0c-4ebf-a8b8-e0f91dff851d
# ╠═54efa70c-bac6-4d7c-93df-0dfd1b89769d
# ╠═5c5d5eeb-6775-452f-880d-7b4fa2acda57
# ╟─61ab2949-d72f-4d80-a717-4b6a9227de0e
# ╠═34b06415-21c1-4904-97f0-ab614447355c
# ╟─a39b90ec-a4c5-472f-b00e-c81ee9c5576f
# ╟─cbbd8670-49ab-4601-b8d7-9f3f456752e8
# ╟─edfbf364-e126-4e95-93d2-a6adfb340045
# ╠═a786ca10-06d2-4b76-97a9-2bcf879ea6cb
# ╟─a3d29aa3-96ca-4681-960c-3b4b04b1e40d
# ╟─6bf06c12-cf25-43c4-81f3-b1d79d13fc94
# ╠═1d72b291-24b7-4ec6-8307-1da0bc4a9183
# ╟─1867ecec-3c3d-4b2b-9036-1488e2184c40
# ╟─0e09d312-2ddf-4d1f-8ad5-a50fb48ca4dd
# ╠═f0ee3633-1d28-4344-9b85-f5433679582a
# ╟─379373b1-e563-4341-9978-5b35c768b5c7
# ╟─8b2725a2-8007-46b8-a160-75042562794d
# ╟─10fa507e-1429-4eb0-b74c-1e6638725690
# ╠═5a3500c2-4f82-43e9-a31b-d530f56fdbe9
# ╟─a1e1d5f8-e06e-4682-ab54-a9454a8e3b30
# ╟─e36979d5-c1b6-4c17-a65a-d8de8e6bd8d0
# ╟─aeafe1ed-f217-48fd-9624-add5f6f791e6
# ╟─b547311c-ddf0-4053-9de4-f0e85b861e63
# ╟─2ed2f3c3-619b-4aed-b88b-b92a43578d84
# ╟─849f69b0-07af-40ab-8295-c0b80a26a2d5
# ╟─36372d31-215d-4299-b4f1-49e42d8b0dbd
# ╟─c32725a4-e276-4372-8d06-d40ba52c9f09
# ╟─f1a71f47-fb19-4988-a439-2ff8d38be5b7
# ╠═a6c4e663-f1e3-4e0c-a8bf-7c13fcb076f0
# ╠═f9e5092b-a158-48cb-a919-96f30ec51038
# ╠═c8c3fe32-560d-11ec-0617-2dc33608384a
# ╠═5394b13a-e629-47c8-902d-685c061b37ae
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
