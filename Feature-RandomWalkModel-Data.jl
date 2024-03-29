### A Pluto.jl notebook ###
# v0.17.4

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

# ╔═╡ 2bb52ee4-1c6f-46b6-b105-86827ada0f75
md"""
# Random Walks and The Efficient Market Hypothesis 
"""

# ╔═╡ 2f499c95-38cf-4856-b199-6c9aac44237a
md"""
### Introduction

[A Random Walk Down Wall Street](https://en.wikipedia.org/wiki/A_Random_Walk_Down_Wall_Street) written by [Princeton University Professor Burton Malkiel](https://dof.princeton.edu/about/clerk-faculty/emeritus/burton-gordon-malkiel) popularized the [efficient market hypothesis (EMH)](https://en.wikipedia.org/wiki/Efficient-market_hypothesis), which posits that stock prices fully reflect all available information and expectations. The efficient market hypothesis was developed by [University of Chicago Professor Eugene Fama](https://en.wikipedia.org/wiki/Eugene_Fama):

* [Fama, E. F. (1970). Efficient Capital Markets: A Review of Theory and Empirical Work. The Journal of Finance, 25(2), 383–417. 	https://doi.org/10.2307/2325486](https://www.jstor.org/stable/2325486)

who was later awarded the [Noble Prize in Economics](https://www.nobelprize.org/prizes/economic-sciences/2013/press-release/) along with  [Robert J. Shiller](http://www.econ.yale.edu/~shiller/) from Yale and [Lars Peter Hansen](https://en.wikipedia.org/wiki/Lars_Peter_Hansen) also from the University of Chicago, for their empirical analysis of asset prices.

Our interest in the efficient market hypothesis follows from its relationship with [random walks](https://en.wikipedia.org/wiki/Random_walk), and in particular, the potential of predicting the price of a stock 𝒯 days into the future (at least in a probabilistic sense) using a random walk approach, an idea first applied to stock prices by [Louis Bachelier](https://en.wikipedia.org/wiki/Louis_Bachelier) in 1900. Fluctuations in prices are explained by the changes in the instantaneous demand and supply of any given stock, causing a random walk in prices. Let's explore random walk models in more detail.
"""

# ╔═╡ 89b5c4d0-68cb-499f-bf99-216e38b40ca0


# ╔═╡ 489464b4-3139-466c-80c7-84449dcec698


# ╔═╡ b4f60458-0299-4246-b2ad-dc3f1db6382b


# ╔═╡ 34bf07e8-5c47-4aa8-aa2c-161709c158be
md"""
### Materials and Methods

##### Disclaimer

This notebook (and all codes discussed within) is offered solely for training and  informational purposes. No offer or solicitation to buy or sell securities or securities derivative products of any kind, or any type of investment or trading advice or strategy,  is made, given, or in any manner endorsed by Pooksoft. 

Trading involves risk. Carefully review your financial situation before investing in securities, futures contracts, options, or commodity interests. Past performance, whether actual or indicated by historical tests of strategies, is no guarantee of future performance or success. Trading is generally not appropriate for someone with limited resources, investment or trading experience, or a low-risk tolerance.  Only risk capital that will not be needed for living expenses.

You are fully responsible for any investment or trading decisions you make, and such decisions should be based solely on your evaluation of your financial circumstances, investment or trading objectives, risk tolerance, and liquidity needs.
"""

# ╔═╡ 65a0683e-3124-4844-ad0a-cccdc9192d05


# ╔═╡ e603c785-b697-48e3-92be-e9213e72215b


# ╔═╡ badff812-6b68-4f6d-b495-3f344873d45c


# ╔═╡ 6a4164c2-bb11-437a-bc7d-a12e363d3e84


# ╔═╡ c81486cb-f78e-4eec-a1cd-1ea113428bd6


# ╔═╡ 880b1174-0925-4e8c-b9f6-f6e295412824
md"""
##### Random Walk Models (RWMs)
[Random walk models](https://en.wikipedia.org/wiki/Random_walk) are tools for computing the time evolution of the prices of risky assets, for example, the price of a stock with the ticker symbol `XYZ`. Let the price of `XYZ` at time $j$ be given by $P_{j}$ (units: USD/share). Then during the next time period (index $j+1$) the `XYZ` share price will be given by:

$$P_{j+1} = P_{j}\exp\left(\mu_{j\rightarrow{j+1}}\Delta{t}\right)$$

where $\Delta{t}$ denotes the length of time between periods $j$ and $j+1$ (units: time) and $\mu_{j\rightarrow{j+1}}$ denotes the growth rate (in the financial world called the [return](https://www.investopedia.com/terms/a/absolutereturn.asp)) between time period(s) $j$ and $j+1$ (units: time$^{-1}$). The growth rate $\mu_{j\rightarrow{j+1}}$ is different for different tickers, and is not constant between time period $j\rightarrow{j+1}$. Dividing both sides by $P_{j}$ and taking the [natural log](https://en.wikipedia.org/wiki/Natural_logarithm) gives the expression:

$$p_{j+1} = p_{j} + \mu_{j\rightarrow{j+1}}\Delta{t}$$

where $p_{\star}$ denotes the log of the price at time $\star$ (units: dimensionless). This expression is a basic [random walk model](https://en.wikipedia.org/wiki/Random_walk). The challenge to using this model is how to estimate the return 
$\mu_{j\rightarrow{j+1}}$ because the values for $\mu_{j\rightarrow{j+1}}$ are [random variables](https://en.wikipedia.org/wiki/Stochastic_process) which are [distributed](https://en.wikipedia.org/wiki/Probability_distribution) in some unknown way. However, if we had a value for $\mu_{j\rightarrow{j+1}}$ (or at least a model to estimate a value for it, denoted as $\hat{\mu}_{j\rightarrow{j+1}}$), we could 
simulate the price of `XYZ` during the next time period $j+1$ given knowledge of the price now (time period $j$). For example, given a value for the close price of `XYZ` on Monday, we could estimate a value for the close price of `XYZ` on Tuesday if we had a model for the daily return.
"""

# ╔═╡ 4c2bbb6b-f287-43c4-9e24-f44fb4c87e36


# ╔═╡ 00120ff5-e9ed-4e71-8208-cf8efd2eac6a


# ╔═╡ e3ba3cf8-2f94-4276-98ed-6fedcfaadd43


# ╔═╡ a64f7082-b834-4400-880a-032cb9aafe4c


# ╔═╡ a1b50b95-bf7d-4f2a-817d-9edb3418aeb3
md"""
##### Estimating models for the return from historical data
Suppose values for the return were governed by some [probability distribution](https://en.wikipedia.org/wiki/Probability_distribution) $\mathcal{D}\left(\bar{m},\sigma\right)$ where $\bar{m}$ 
denotes the [mean value](https://en.wikipedia.org/wiki/Mean) and $\sigma$ denotes the [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) of the growth rate $\mu_{j\rightarrow{j+1}}$. Of course, we do not know the actual values for 
$\left(\bar{m},\sigma\right)$, nor do we know the form of $\mathcal{D}\left(\bar{m},\sigma\right)$, but we can estimate them from historial price data. 

An even better way to learn $\mathcal{D}\left(\bar{m},\sigma\right)$ would be to estimate the form of the distribution from the data itself using a technique such as [Kernel Density Estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation) or KDE. 
"""

# ╔═╡ 1852d4c5-e73d-4038-a565-b8fb6ff63502


# ╔═╡ 1c9f096f-7dce-4d0c-8a71-6fbe682e514d


# ╔═╡ eb4a5874-c13c-4915-96b9-0c47eaa7f50c


# ╔═╡ 66bfd748-3b40-42df-86fa-e55019dba856


# ╔═╡ cf11b553-13e5-488d-bf1c-16de5911d658


# ╔═╡ 2353a285-71de-43b2-a60f-5a3274ff9e6b
md"""
##### My Julia Setup
In the following code block, we set up project paths and import external packages that we use to download and analyze historical price history for a variety of ticker symbols. In addition, we load a local copy of the [Serenity library](https://github.com/Pooksoft/Serenity).
"""

# ╔═╡ 43d75d79-b710-4dc5-9478-dd2b08616be9


# ╔═╡ ebc5ed32-1fe3-4854-a326-7a068e14164b


# ╔═╡ 1b25e7d1-909f-4fdc-9700-5d131251e1b5


# ╔═╡ 9943d000-83d0-413d-a231-0295fb19df71
md"""
### Results and Discussion
"""

# ╔═╡ f66a480b-3f0c-4ebf-a8b8-e0f91dff851d
md"""
##### Download historical price data from the [Polygon.io](https://www.polygon.io) financial data warehouse
In this study, we model two years of daily adjusted close price data from `2019-12-26` to `2021-12-23` (N = 504) for a variety of ticker symbols (𝒫 = 40). We used methods in the [Serenity library](https://github.com/Pooksoft/Serenity) and the routine `aggregates_api_endpoint` defined in this notebook to download historical pricing data from the [Polygon.io](https://www.polygon.io) financial data warehouse using a free tier application programming interface key.

The `aggregates_api_endpoint` takes a list of ticker symbols that we want to model as an input argument.
This method checks to see if we have price data already saved locally for each ticker.

* If yes, then we load the saved file as a [DataFrame](https://dataframes.juliadata.org/stable/) and store it in the `price_data_dictionary` data structure which is of type [Dict](https://docs.julialang.org/en/v1/base/collections/#Dictionaries) where the keys are the ticker strings e.g., MSFT, etc and the values are the price [DataFrames](https://dataframes.juliadata.org/stable/).
* if no, we download new data from [Polygon.io](https://www.polygon.io), save this data locally as a `CSV` file (<ticker>.csv), and return the price data as a [DataFrame](https://dataframes.juliadata.org/stable/) in the `price_data_dictionary` dictionary.

The parameters for the [Polygon.io](https://www.polygon.io) application programming interface call are encoded in endpoint-specific [Structs](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) defined in the [Serenity library](https://github.com/Pooksoft/Serenity). 
"""

# ╔═╡ f94ecbae-e823-46f6-a69f-eb7edce7cfe5


# ╔═╡ dc5f6727-cd17-4c75-80ce-94915a0e359a


# ╔═╡ 9b9d7fdf-91b1-46e4-bd66-43e50add56be


# ╔═╡ 0e7d1741-88a1-4e8e-b964-e8ead4d1807e


# ╔═╡ 300b62c8-830a-472f-a07c-17153468c1fb
md"""
##### What ticker symbols are we going to model?
"""

# ╔═╡ 54efa70c-bac6-4d7c-93df-0dfd1b89769d
# Pooksoft Industrial Average (PSIA) -> the DJIA + some stuff
ticker_symbol_array = sort([
	"MSFT", "ALLY", "MET", "AAPL", "GM", "PFE", "TGT", "WFC", "AIG", "F", "GE", "AMD",
	"MMM", "AXP", "AMGN", "BA", "CAT", "CVX", "CSCO", "KO", "DIS", "DOW", "GS", "HD", "IBM",
	"HON","INTC", "JNJ", "JPM", "MCD", "MRK", "NKE", "PG", "CRM", "TRV", "UNH", "VZ", 
	"V", "WBA", "WMT", "SPY", "SPYD", "SPYG", "SPYV", "SPYX", "VOO", "VTI", "VEA", "VWO",
	"VNQ", "VGK", "MRNA"
]);

# ╔═╡ 866cc84d-86c1-40e2-bd29-deae01da9a2e
𝒫 = length(ticker_symbol_array); # the number of ticker symbols is given the symbol 𝒫

# ╔═╡ 7bcbdc4e-a38a-4201-a1ec-d2e4df4d2f6a
𝒯 = 14; # number of days that we simulate in the future

# ╔═╡ d1edad45-5df3-43cd-8abc-97d84fab699b


# ╔═╡ bb7582b4-3ecc-44a7-810a-aae0dc2fe816


# ╔═╡ 34db9c1e-6af7-4710-8de4-fd0caacbde36


# ╔═╡ 46fef026-6bac-401b-ab6f-75f2aff79e6a


# ╔═╡ 455b2bea-4d79-4dd8-964c-80835bb88727


# ╔═╡ 4d4230e5-28dc-4f87-b732-686a5c91fc4f
md"""
##### Pull down data from Polygon.io for our list of ticker symbols
"""

# ╔═╡ e2eb17d0-2ca3-437c-97f7-04e76ee879cb


# ╔═╡ a3a08666-f494-489d-a888-8fe7baa70729


# ╔═╡ b05cac4c-d5d9-47f7-ab38-d9be8322ab45


# ╔═╡ 61ab2949-d72f-4d80-a717-4b6a9227de0e
md"""
##### Estimate the daily return distributions

To estimate the return moel $\mathcal{D}\left(\bar{m},\sigma\right)$, we compute the historical returns for each ticker using the `compute_log_return_array` method in the [Serenity library](https://github.com/Pooksoft/Serenity.git).

Given the returns, we estimate the parameters of $\mathcal{D}\left(\bar{m},\sigma\right)$ using a variety of approaches e.g., [nonlinear least squares](https://en.wikipedia.org/wiki/Non-linear_least_squares) or [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation). However, there are a few technical questions with this approach:

* __Question 1__: What model for $\mathcal{D}\left(\bar{m},\sigma\right)$ do we choose? Conventional wisdom suggests a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution). However, is this really correct, or is the data governed by another type of distribution?
* __Question 2__: An underlying assumption of the random walk is that the daily returns are independent and identically distributed, ie., the close price for the days $j$ and $j+1$ are independent of one another, they follow the same underlying distribution, and the parameters $(\bar{m},\sigma)$ don't change in time.
"""

# ╔═╡ e461b560-435b-4104-b42e-af04b1d25984


# ╔═╡ b7b618c2-d48b-45c7-aff8-600eda010169


# ╔═╡ a39b90ec-a4c5-472f-b00e-c81ee9c5576f
md"""
__Fig. 1__: Histogram of the adjusted daily returns computed using the `stephist` routine of the [StatsPlots.jl package](https://github.com/JuliaPlots/StatsPlots.jl) for 𝒫 = 40 ticker symbols in PSIA. The number of bins was set to 20% of the number of historical records for `AAPL` (shown in red) and was constant for each ticker. The returns were calculated using approximately two years of historical daily adjusted close price data (depending upon the ticker symbol). 
Close price data was downloaded from the [Polygon.io](https://www.polygon.io) data warehouse. 
"""

# ╔═╡ 3f9879fa-ec04-4d13-8207-2c77d9ddfca2


# ╔═╡ 22992de2-dfba-4a18-8909-62bddbf4e0e7


# ╔═╡ 367158be-9a7f-4b76-96c9-2806f9aad75e


# ╔═╡ fb83174c-fcde-4496-8232-545f17ac9d2d


# ╔═╡ edfbf364-e126-4e95-93d2-a6adfb340045
md"""
##### The daily return of individual ticker symbols is not normally distributed

A key assumption often made in mathematical finance is that the return $\mu_{j\rightarrow{j+1}}$ (growth rate) computed using the log of the prices is normally distributed (Question 1 raised above). However, a closer examination of the return data for our collection of ticker symbols suggests this is not true.

The return $\mu_{j\rightarrow{j+1}}$ values computed from historical data are not normally distributed for the vast majority of ticker symbols in our collection (Fig. 2). Visual inspection of the histograms computed using a [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution) (Fig 2, dark blue) for the returns appear to be closer to the historical data (Fig 2, red) compared to a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) (Fig 2, light blue). [QQplots](https://en.wikipedia.org/wiki/Q–Q_plot) and the [Kolmogorov–Smirnov (KS) test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test) further support this finding. 
"""

# ╔═╡ a786ca10-06d2-4b76-97a9-2bcf879ea6cb
# fit a distribution to a ticker -
single_asset_ticker_symbol = "GS";

# ╔═╡ f067b633-fde0-4743-bd76-f8c390a90950


# ╔═╡ a3d29aa3-96ca-4681-960c-3b4b04b1e40d
md"""
__Fig 2__: Comparison of the actual (red line) and estimated histogram for the adjusted daily returns for ticker = $(single_asset_ticker_symbol) using a [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution) for the return model (blue line) and a [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) for the return model (light blue line). The return model distributions were estimated using the `fit` routine of the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package. The `fit` routine uses multiple approaches to estimate the parameters in the distribution including [Maximum Likelihood Estimation (MLE)](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation). 
"""

# ╔═╡ 4dff7597-bc0e-4288-9413-c26a132c1e44


# ╔═╡ 1d72b291-24b7-4ec6-8307-1da0bc4a9183
md"""
Beyond simply visualizing the histogram, we explored the question of a [Normal](https://en.wikipedia.org/wiki/Normal_distribution) versus [Laplace](https://en.wikipedia.org/wiki/Laplace_distribution) return distribution, by constructing [QQplots](https://en.wikipedia.org/wiki/Q–Q_plot) to visualize the historical return data versus a Laplace distribution (Fig. 3). Further, we performed the [Kolmogorov–Smirnov (KS) test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test), using the KS test implementation in the [HypothesisTests.jl](https://juliastats.org/HypothesisTests.jl/latest/nonparametric/#Kolmogorov-Smirnov-test-1) package. The KS test determines whether data follows a particular distribution, in our case we test whether the historical return data were either Normally or Laplace distributed (Table 1). 
"""

# ╔═╡ 1867ecec-3c3d-4b2b-9036-1488e2184c40
md"""
__Fig 3__: Quantile-Quantile plot (QQplot) for historical $(single_asset_ticker_symbol) return data versus theoretical data generated using a Laplace distribution. If the historical data were governed by a Laplace distribution, the points would lie along the equality line. 
"""

# ╔═╡ a3f1710e-eb98-46c6-aadb-5cf0b98e1bc6


# ╔═╡ 9380908a-8cc5-4d3a-9d6c-03b491198bc1


# ╔═╡ a07d661a-a0b4-4000-b7a4-5f17cda5edc3
md"""
###### Example:  Kolmogorov–Smirnov (KS) test for a Laplace and Normal distrubution for ticker: $(single_asset_ticker_symbol)

The null hypothesis $h_0$ for the KS test is the data follows a specified distribution, in our case a Laplace distribution. Conversely, the alternative hypothesis $h_{1}$ is the data does not follow a particular distribution, again in this case, a Laplace distribution. 
"""

# ╔═╡ fdc34643-2c27-4bc9-a077-8bac68dc56b7


# ╔═╡ 8897de5b-a6c7-4b05-98c6-d738bbb527af


# ╔═╡ 379373b1-e563-4341-9978-5b35c768b5c7
md"""
__Table 1__: [Kolmogorov–Smirnov (KS)](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test) test results for Normal and Laplace distribution for PSIA tickers (𝒫 = 40). The historical return values computed over the last two years of historical data for 39 of the 40 members of the PSIA were governed by a Laplace distribution; MRNA failed the KS test for a Laplace distribution, and was estimated to be Normally distributed.  
"""

# ╔═╡ 95495255-385a-424b-9b84-dcd8a1187282


# ╔═╡ d8af95e4-fe4c-4cf0-8195-dab80629177c


# ╔═╡ 3f62e075-4724-456e-ab59-9b85d4263ee6


# ╔═╡ bc4e93f5-7e76-4115-8cfd-039212181fb3
md"""
##### Are daily returns iid? The Wald-Wolfowitz runs test

A second key assumption of the random walk approach is that returns are independently and identically distributed. To test whether this is true for the ticker symbols in our data set, we performed the [Wald-Wolfowitz runs test](https://en.wikipedia.org/wiki/Wald–Wolfowitz_runs_test). The [Wald-Wolfowitz runs test](https://en.wikipedia.org/wiki/Wald–Wolfowitz_runs_test) is a non-parametric tool, that does not require returns to be normally distributed, which calculates the likelihood that a sequence is independent (the null hypothesis). 

Each return is scored according to whether it is larger, smaller, or the same as the median; a score of +1 is assigned when the return is greater than the median, a score of -1 is assigned when the return is less than the median, and a score of 0 is assigned when the return equals the mean. 
"""

# ╔═╡ 587c70c4-a3b8-4300-b7a4-aa9f6b1bc276


# ╔═╡ 421e213e-9780-4a4d-a411-009541c44a9e


# ╔═╡ 2d40603b-56e3-49ca-a8ea-e13c933f5e19
md"""
__Table XX__: Wald-Wolfowitz test results for historical price data from `2019-12-26` to `2021-12-23` for 𝒫 = 40 tickers in the PSIA. The Wald-Wolfowitz test implementation from the [HypothesisTests.jl]() package was used to compute the test results.
"""

# ╔═╡ c39168a1-dee2-4421-8f44-266df47dd08e


# ╔═╡ 4e6e394d-0d95-46cc-8344-0222b0fb761e


# ╔═╡ 0454a796-b35d-4caa-b847-f578191f2896


# ╔═╡ 10fa507e-1429-4eb0-b74c-1e6638725690
md"""
##### Monte Carlo simulations of the daily close price using a random walk model

We used [Monte-Carlo simulation](https://www.investopedia.com/articles/investing/112514/monte-carlo-simulation-basics.asp) to estimate the close price distribution of each ticker, 𝒯 days into the future, using the daily return model $\mathcal{D}\left(\bar{m},\sigma\right)$ we learned from historical price data. Monte-Carlo simulation is conceptually easy to understand but tricky to implement in an efficient way. 

Imagine that we split our reality into a large number of structurally similar but parallel dimensions, such that in each of these dimensions there was a stock market in which company `XYZ` traded, and for each dimension, yesterday's close price for `XYZ` was $p_{j}$. Next, further, suppose that we cloned ourselves and placed a clone into each of these dimensions so that we could record the stock price of ticker `XYZ` at the end of the day for the specified number of days. Because dimensions were not identical, we expect the close price to vary between dimensions. After 𝒯 days the clones report their price recordings we put them all together. The data provided by each clone is called a _sample path_, and we can use these sample paths to compute a distribution of possible outcomes for the price of `XYZ`.  

In a technical sense, starting from close price at time index $p_j$, we draw a random sample from the $\mathcal{D}\left(\bar{m},\sigma\right)$ distribution, and compute a value for $p_{j+1}$. We then advance the time index $p_{j+1}\rightarrow{p_{j}}$ and repeat the process over again until we reach the end of the prediction time horizon. Of course, because we are dealing with random variables, we repeat the calculation of each sample path many times which gives a distribution of possible price outcomes. 

"""

# ╔═╡ 92673e3f-e436-4c7e-accf-216574233d58


# ╔═╡ e2df85a5-975c-472b-8538-d07b682d1647


# ╔═╡ 2e416e6d-ea7c-464c-9632-83b0f7f06b2a


# ╔═╡ 3d1a1f06-5707-46ff-b44a-539c62fe008b


# ╔═╡ 90a0a645-cd70-4edb-8240-a152d6e7bb3a


# ╔═╡ cdd66194-cb79-433e-a16c-64be76de83a4


# ╔═╡ 4f8a3476-93a8-414e-b169-7046f1a57547


# ╔═╡ 21791e69-db82-4081-97ec-e7f2c0a46b5d


# ╔═╡ 1f9459dd-57ba-4353-81b2-9e04aa8257af


# ╔═╡ 90722894-dc53-4964-8141-e35adfeb8a76


# ╔═╡ 4d97d757-ed9f-4e77-9b54-6b9c8f2f49ee


# ╔═╡ 9a3d5138-0490-42db-9f99-310d94124951


# ╔═╡ bdf1e4d0-d2e6-4ab5-92c5-5f25518a9acc


# ╔═╡ 0a55c3a1-a834-4ec4-beac-0788091f7a70


# ╔═╡ 4725bf3c-f058-4238-8c8f-297dc2851c6e


# ╔═╡ 3b9d8429-375d-46c3-a115-eb9022262e09


# ╔═╡ 23ae5b4c-f690-46a6-b2c9-3d753cf247d6


# ╔═╡ edf426b6-a571-4938-890a-01a089d02b29


# ╔═╡ c32725a4-e276-4372-8d06-d40ba52c9f09
md"""
### Conclusions

In this study, we explored random walks and the efficient market hypothesis. In particular, we simulated 𝒯 future price predictions based on a return model learned from historical data. There were a few takeaway points:

* Daily returns were not Normally distributed for a variety of ticker symbols (the DJIA plus 10 additional symbols) with the exception of MRNA, which appears to be normally distributed. Ticker returns more likely followed a Laplace distribution.
* The Wald-Wolfowitz runs test suggested that the daily returns for several members of the DJIA were not random, thus, violating the iid assumption of our random walk modeling approach. However, this violation did not appear to hinder the Monte-Carlo estimate of the close price distribution. 
* A random walk model predicted future stock price distributions for short time horizons (𝒯 = 7 days), but was less successful for longer time horizons e.g., 𝒯 = 35 days.

"""

# ╔═╡ f4bd98b5-f4a8-424a-b974-41aceede92fb


# ╔═╡ 6086ed53-4fbe-4037-b435-8d8aa20a1417


# ╔═╡ defa9be8-4d19-488b-820c-2a3526401cbf


# ╔═╡ f1a71f47-fb19-4988-a439-2ff8d38be5b7
function ingredients(path::String)
	
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol("Serenity")
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ fe2848df-823a-4ed0-918c-2c200957ee80
begin
	
	# setup: we need to specify the paths where we can find our project resources
    _PATH_TO_ROOT = pwd()
    _PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")
    _PATH_TO_DATA = joinpath(_PATH_TO_ROOT, "data")
    _PATH_TO_CONFIG = joinpath(_PATH_TO_ROOT, "configuration")

    # what packages are we going to use?
    using TOML
	using PlutoUI
	using HypothesisTests
	using PrettyTables
	using Colors
	using StatsPlots
	using StatsBase
	using Reexport
	
	# these packages are reexported by Serenity 
	using DataFrames
    using CSV
    using HTTP
    using Dates
    using Statistics
    using Distributions
	using Optim
	using Convex
	using SCS
	using MathOptInterface
	using JSON
	using Colors
	
    # alias the Polygon.io URL -
    DATASTORE_URL_STRING = "https://api.polygon.io"

    # What we have here is a classic good news, bad news situation ...
    # Bad news: We don't check in our Polygon.io API key to GitHub (sorry). 
    # Good news: Polygon.io API keys are free. 
    # check out: https://polygon.io -
    configuration_dictionary = TOML.parsefile(joinpath(_PATH_TO_CONFIG, "Configuration.toml"))

    # get some stuff from the configuration dictionary -
    POLYGON_API_KEY = configuration_dictionary["API"]["polygon_api_key"]

	# load the local version of the serenity library -
	Serenity = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"));

    # show -
    nothing
end

# ╔═╡ e36979d5-c1b6-4c17-a65a-d8de8e6bd8d0
md"""
Show actual price trajectory? $(@bind show_real_traj CheckBox()) 
"""

# ╔═╡ b04cba56-dd48-403b-82fc-1cf3713853a7
begin

	# List of colors -
	WHITE = RGB(1.0, 1.0, 1.0)
	BACKGROUND = RGB(0.99, 0.98, 0.96)
	BLUE = RGB(68 / 255, 119 / 255, 170 / 255)
	LBLUE = RGB(102 / 255, 204 / 255, 238 / 255)
	GRAY = RGB(187 / 255, 187 / 255, 187 / 255)
	RED = RGB(238 / 255, 102 / 255, 119 / 255)

	# show -
	nothing
end

# ╔═╡ 91dae79f-e454-4b27-84a7-4cbc6bc33265
function aggregates_api_endpoint(ticker_symbol::String; sleeptime::Float64 = 20.0)::NamedTuple

	# save the file locally to avoid the API limit -
	savepath = joinpath(_PATH_TO_DATA,"polygon", "random_walk_model", "aggregates")
	local_path_to_data_file = joinpath(savepath, "$(ticker_symbol).csv")
	local_path_to_header_file = joinpath(savepath, "header", "$(ticker_symbol)-header.csv")

	# do we have a saved file already?
	if (ispath(local_path_to_data_file) == true && 	ispath(local_path_to_header_file) == true)
	
		# we already have this ticker downloaded. load the existing file from disk -
        data_table = CSV.read(local_path_to_data_file, DataFrame)
		header = CSV.read(local_path_to_header_file, Dict)

        # put the price DataFrame into the price_data_dictionary -
        return (header=header, data=data_table)
		
	else

		# make the dir to save the data (if we don't have it already) -
		mkpath(savepath)
    	mkpath(joinpath(savepath,"header"))
		
		# Build an API model for the aggregates endpoint -
		# see: https://polygon.io/docs/stocks/getting-started
		aggregates_api_model = Serenity.PolygonAggregatesEndpointModel()
		aggregates_api_model.adjusted = true
		aggregates_api_model.sortdirection = "asc"
		aggregates_api_model.apikey = POLYGON_API_KEY
		aggregates_api_model.limit = 5000
		aggregates_api_model.to = Date(2021, 12, 24)
		aggregates_api_model.from = Date(2019,12, 24)
		aggregates_api_model.multiplier = 1
		aggregates_api_model.timespan = "day"
		aggregates_api_model.ticker = ticker_symbol

		# execute the API call for these parameters -
		(header,df) = Serenity.polygon(DATASTORE_URL_STRING, aggregates_api_model);

		# save locally -
		CSV.write(local_path_to_data_file, df)
		CSV.write(local_path_to_header_file, header)

		# sleep sometime before we return - helps avoid API limit issue -
		sleep(sleeptime)
			
        # put the price DataFrame into the price_data_dictionary -
        return (header=header, data=df)
	end
end

# ╔═╡ 5c5d5eeb-6775-452f-880d-7b4fa2acda57
begin

	# initialize -
	price_data_dictionary = Dict{String,DataFrame}()

	# process each member of the ticker symbol array -
	for ticker_symbol ∈ ticker_symbol_array

		# make a call to Polygon.io (or load local cached data)
		# the aggregates_api_endpoint method is defined at the bottom of this notebook
		(h,d) = aggregates_api_endpoint(ticker_symbol)

		# grab the price data -
		price_data_dictionary[ticker_symbol] = d
	end

	# show -
	nothing
end

# ╔═╡ 34b06415-21c1-4904-97f0-ab614447355c
begin
    # initialize -
    return_data_dictionary = Dict{String,DataFrame}()

    # compute -
    for ticker_symbol ∈ ticker_symbol_array

        # Serenity -> computes the log return given a DataFrame, returns a DataFrame
        return_data_dictionary[ticker_symbol] =
			Serenity.compute_log_return_array(price_data_dictionary[ticker_symbol],
            :timestamp => :close; Δt = (1.0 / 1.0))
    end
end

# ╔═╡ cbbd8670-49ab-4601-b8d7-9f3f456752e8
begin

    # compute some ranges -
    ticker_symbol = "AAPL"
    μ_avg = mean(return_data_dictionary[ticker_symbol][!, :μ])

    # number of bins - 10% of the length of a test ticker
    number_of_bins = convert(Int64, (floor(0.20 * length(return_data_dictionary[ticker_symbol][!, :μ]))))
	
    # number of tickers -
    number_of_ticker_symbols = length(ticker_symbol_array)
    for ticker_index = number_of_ticker_symbols:-1:1
        local_ticker_symbol = ticker_symbol_array[ticker_index]

        if (ticker_index == number_of_ticker_symbols)

            # make the first plot -
            stephist(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true,
                background_color = BACKGROUND, background_color_outside = WHITE,
                foreground_color_minor_grid = RGB(1.0, 1.0, 1.0),
                lw = 1, c = GRAY, foreground_color_legend = nothing, label = "")

        elseif (ticker_index != 1 && ticker_index != number_of_ticker_symbols)

            stephist!(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true, lw = 1,
                c = GRAY, label = "")
        else
            stephist!(return_data_dictionary[local_ticker_symbol][!, :μ], bins = number_of_bins, normed = :true, lw = 2,
                c = RED, label = "$(local_ticker_symbol)")
        end
    end

    # label the plots -
    xlabel!("Daily return μ (1/day)", fontsize = 18)
    ylabel!("Frequency (N=$(number_of_bins); dimensionless)", fontsize = 14)
    xlims!((-100.0 * μ_avg, 100.0 * μ_avg))
end

# ╔═╡ 012ef00f-6176-4d42-801f-5765c47df7ba
begin 

	# get the historical return data -
    μ_vector = return_data_dictionary[single_asset_ticker_symbol][!, :μ]

	# fit a Laplace -
	LAPLACE = fit(Laplace, μ_vector)

	# fit a normal -
    NORMAL = fit(Normal, μ_vector)

	# show -
	nothing
end

# ╔═╡ f0ee3633-1d28-4344-9b85-f5433679582a
let
	# perform the KS test -
	test_result = HypothesisTests.ExactOneSampleKSTest(unique(μ_vector), LAPLACE)
end

# ╔═╡ 345d8feb-cfd0-4acb-a626-dc16181ddb68
let
	# perform the KS test for Laplace -
	test_result = HypothesisTests.ExactOneSampleKSTest(unique(μ_vector), NORMAL)
end

# ╔═╡ 6bf06c12-cf25-43c4-81f3-b1d79d13fc94
let
    
    # generate 20_000 samples
    S = rand(LAPLACE, 20000)
    SN = rand(NORMAL, 20000)

    # plot against actual -
    stephist(return_data_dictionary[single_asset_ticker_symbol][!, :μ],
		bins = convert(Int64,round(1.2*number_of_bins)), 
		normed = :true, 
		lw = 2, c = RED,
        label = "$(single_asset_ticker_symbol)", background_color = BACKGROUND,
        background_color_outside = WHITE, foreground_color_legend = nothing)

    stephist!(S, bins = number_of_bins, normed = :true, lw = 2, c = BLUE,
        label = "$(single_asset_ticker_symbol) Laplace Model")

    stephist!(SN, bins = number_of_bins, normed = :true, lw = 2, c = LBLUE,
        label = "$(single_asset_ticker_symbol) Normal Model")

    xlabel!("Daily return μ (1/day)", fontsize = 18)
    ylabel!("Frequency (N=$(number_of_bins); dimensionless)", fontsize = 14)
    xlims!((-75.0 * μ_avg, 75.0 * μ_avg))
	
end

# ╔═╡ 0e09d312-2ddf-4d1f-8ad5-a50fb48ca4dd
let

	# get the historical return data -
    μ_vector = return_data_dictionary[single_asset_ticker_symbol][!, :μ]
	
	qqplot(Laplace, μ_vector, label="QQPlot $(single_asset_ticker_symbol)", legend=:topleft, 
		foreground_color_legend = nothing, c=LBLUE, markerstrokecolor=BLUE, lw=2, 
		background_color = BACKGROUND, background_color_outside = WHITE)
	xlabel!("Theoretical Laplace Quantiles", fontsize=18)
	ylabel!("Sample Quantiles", fontsize=18)
end

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

# ╔═╡ 20bba789-504b-4d27-a2fc-4f08badc9a53
begin
	μ_local = return_data_dictionary[single_asset_ticker_symbol][!,:μ]
	WaldWolfowitzTest(μ_local)
end

# ╔═╡ 1c816d04-3b76-4c0c-9ed9-9b20b5ad50d9
with_terminal() do

	# initialize -
	run_array = Array{Union{Int,String, Float64},2}(undef,𝒫,10)
	
	for (ticker_index, ticker) ∈ enumerate(ticker_symbol_array)

		# get the return for this ticker -
		μ_local = return_data_dictionary[ticker][!,:μ]
		N = length(μ_local)

		# what is the mean return for this ticker?
		median_value = median(μ_local)
		tmp_array = Array{Int64,2}(undef,N,3)
		fill!(tmp_array,0)
		for time_index ∈ 1:N
		
			# classify -
			if (μ_local[time_index] >= median_value)
				tmp_array[time_index,1] += 1
				tmp_array[time_index,2] += 0
				tmp_array[time_index,3] += 1
			elseif (μ_local[time_index] < median_value)
				tmp_array[time_index,1] += 0
				tmp_array[time_index,2] += 1
				tmp_array[time_index,3] = -1
			end
		end

		run_array[ticker_index,1] = ticker
		run_array[ticker_index,2] = N

		# count the number of runs -
		U = 1
		for time_index ∈ 2:N
			old_value = tmp_array[time_index - 1,3]
			new_value = tmp_array[time_index,3]
			if (old_value!=new_value)
				U += 1
			end
		end

		run_array[ticker_index,3] = sum(tmp_array[:,1])
		run_array[ticker_index,4] = sum(tmp_array[:,2])
		run_array[ticker_index,5] = U

		# run the test -
		local_test_result = WaldWolfowitzTest(μ_local)
		run_array[ticker_index,6] = local_test_result.μ
		run_array[ticker_index,7] = local_test_result.σ
		run_array[ticker_index,8] = local_test_result.z
		run_array[ticker_index,9] = pvalue(local_test_result)
		run_array[ticker_index,10] = (pvalue(local_test_result) >= 0.05) ? 1 : 0
	end

	header_data = (["ticker", "N", "N₊", "N₋", "U","μ","σ","Z","p-value", "p-value ≥ 0.05"])
	pretty_table(run_array, header=header_data)
end

# ╔═╡ 5a3500c2-4f82-43e9-a31b-d530f56fdbe9
begin
    
	# how many steps, sample paths etc -
    number_of_sample_paths = 12500;
	number_of_strata = 1;
	monte_carlo_simulation_dictionary = Dict{String,Array{Float64,2}}()

	for ticker_symbol ∈ ticker_symbol_array

		# what is the *actual* price data?
    	local actual_price_data = price_data_dictionary[ticker_symbol][end-𝒯:end, :close]

    	# get initial price -
    	initial_price_value = log(actual_price_data[1])

    	# compute a set of possible trajectories -> convert back to actual price -
    	monte_carlo_simulation_dictionary[ticker_symbol] = LAPLACE(initial_price_value, 𝒯;
        		number_of_sample_paths = number_of_sample_paths, number_of_strata = number_of_strata) .|> exp
	end
    
    # show -
    nothing
end

# ╔═╡ ff11bdd2-8ab6-40c0-844b-a3474d7e9a04
Z = monte_carlo_simulation_dictionary[single_asset_ticker_symbol]

# ╔═╡ cac02388-31c4-40b9-8288-a6685e1854fb
begin
	
	# plot the histogram -
	stephist(Z[end,:], normed=true, lw=2, c=BLUE, 
		label="Close price $(single_asset_ticker_symbol) (T = $(𝒯))", 
		background_color = BACKGROUND, background_color_outside = WHITE, 
		foreground_color_legend = nothing)
	xlabel!("Simulated close price $(single_asset_ticker_symbol) (T = $(𝒯)) (USD/share)", fontsize=18)
	ylabel!("Frequency (dimensionless)", fontsize = 14)	
end

# ╔═╡ a1e1d5f8-e06e-4682-ab54-a9454a8e3b30
md"""
__Fig 4__: In sample random walk simulation of ticker = $(single_asset_ticker_symbol) for a 𝒯 = $(𝒯) day prediction horizon. Blue lines denotes simulated sample paths while the red line denotes the actual price trajectory for ticker $(single_asset_ticker_symbol). The simulation consisted of N = $(2*number_of_sample_paths) sample paths.
"""

# ╔═╡ b547311c-ddf0-4053-9de4-f0e85b861e63
md"""
__Table 2__: Comparison of the actual versus simulated close price for a 𝒯 = $(𝒯) day prediction horizon for each ticker in the PSIA (𝒫 = 40). Each ticker was classified c ∈ {-1,0,1} based upon whether the actual close price Pₐ ∈ Pₑ ± σ, where Pₐ denotes the actual close price (units: USD/share), Pₑ denotes the mean simulated close price (units: USD/share), and σ denotes the standard deviation of the simulated close price (units: USD/share) computed over the family Monte Carlo trajectories (N = $(2*number_of_sample_paths)). Classes: +1 HIGH, 0 INSIDE, or -1 LOW. 
"""

# ╔═╡ aeafe1ed-f217-48fd-9624-add5f6f791e6
begin
	
	# setup some preliminaries -
	skip_factor = convert(Int64, round(0.01*(2*number_of_sample_paths*number_of_strata)))
    plot_index_array = 1:skip_factor:(2*number_of_sample_paths*number_of_strata) |> collect
	simulated_price_trajectory = monte_carlo_simulation_dictionary[single_asset_ticker_symbol]
	
	# what is the *actual* price data?
    actual_price_data = price_data_dictionary[single_asset_ticker_symbol][end-𝒯:end, :close]
	
    # plot -
    plot(simulated_price_trajectory[:, plot_index_array], c = LBLUE, legend = false, label = "", lw = 1,
        background_color = BACKGROUND, background_color_outside = WHITE)

    if (show_real_traj == true)
        plot!(actual_price_data[1:end], c = RED, lw = 3, legend = :topleft, 
			label = "$(single_asset_ticker_symbol) actual", foreground_color_legend = nothing)
    end

    xlabel!("Time step index (day)", fontsize = 18)
    ylabel!("Simulated $(single_asset_ticker_symbol) close price (USD/share)", fontsize = 14)
    # title!("Random walk simulation $(single_asset_ticker_symbol) (N = $(number_of_sample_paths))", fontsize=12)
end

# ╔═╡ 2ed2f3c3-619b-4aed-b88b-b92a43578d84
with_terminal() do

	# initialize some storage -
	price_state_table = Array{Any,2}(undef, number_of_ticker_symbols, 10)
	
	for ticker_symbol_index ∈ 1:number_of_ticker_symbols 

		# get the symbol -
		ticker_symbol = ticker_symbol_array[ticker_symbol_index]

		# compute some data re this symbol -
		simulated_price_trajectory = monte_carlo_simulation_dictionary[ticker_symbol]
		estimated_mean_price = round(mean(simulated_price_trajectory[end, :]), sigdigits = 4)
    	std_estimated_price = round(std(simulated_price_trajectory[end, :]), sigdigits = 4)
		actual_price_data = price_data_dictionary[ticker_symbol][end-𝒯:end, :close]
    	price_actual = actual_price_data[end]
		tmp_LB = estimated_mean_price - 1*std_estimated_price
        tmp_UB = estimated_mean_price + 1*std_estimated_price
		
		# populate state table -
		price_state_table[ticker_symbol_index,1] = ticker_symbol
		price_state_table[ticker_symbol_index,2] = actual_price_data[1]
		price_state_table[ticker_symbol_index,3] = price_actual
		price_state_table[ticker_symbol_index,4] = estimated_mean_price
		price_state_table[ticker_symbol_index,5] = ((price_actual - estimated_mean_price)/price_actual)*100
		price_state_table[ticker_symbol_index,6] = tmp_LB
		price_state_table[ticker_symbol_index,7] = tmp_UB
		
        δL = max(0,round(tmp_LB - price_actual, sigdigits = 4))
        δU = max(0,round(price_actual - tmp_UB, sigdigits = 4))
		price_state_table[ticker_symbol_index,8] = δL
		price_state_table[ticker_symbol_index,9] = δU	

		if (δL == 0.0 && δU == 0.0)
			price_state_table[ticker_symbol_index,10] = 0 	# exepcted
		elseif (δL>0.0 && δU == 0.0)
			price_state_table[ticker_symbol_index,10] = -1 # oversold
		elseif (δL == 0.0 && δU > 0.0)
			price_state_table[ticker_symbol_index,10] = 1 # overbought
		end
	end


	price_table_header = (
		["ticker", "P₀ (𝒯 = 0)", "Pₐ (𝒯 = 14)", "Pₑ (𝒯 = 14)", "ΔP/Pₐ", "ℒ = (Pₑ - σ)", "𝒰 = (Pₑ + σ)", "δL = max(0, ℒ - Pₐ)", 
			"δU = max(0, Pₐ - 𝒰)", "class c"],
		["", "USD/share", "USD/share", "USD/share", "%", "USD/share", "", "USD/share", "USD/share","c ∈ {-1,0,1}"]
	)

	pretty_table(price_state_table, header=price_table_header)
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
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Convex = "f65535da-76fb-5f13-bab9-19810c17039a"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Reexport = "189a3867-3050-52da-a836-e630ba90ab69"
SCS = "c946c3f1-0d1f-5ce8-9dea-7daa1f7e2d13"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[compat]
CSV = "~0.9.11"
Colors = "~0.12.8"
Convex = "~0.14.18"
DataFrames = "~1.3.1"
Distributions = "~0.25.37"
HTTP = "~0.9.17"
HypothesisTests = "~0.10.6"
JSON = "~0.21.2"
MathOptInterface = "~0.10.6"
Optim = "~1.6.0"
PlutoUI = "~0.7.27"
PrettyTables = "~1.3.1"
Reexport = "~1.2.2"
SCS = "~0.8.1"
StatsBase = "~0.33.13"
StatsPlots = "~0.14.30"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AMD]]
deps = ["Libdl", "LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "fc66ffc5cff568936649445f58a55b81eaf9592c"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.4.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

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

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "940001114a0147b6e4d10624276d56d531dd9b49"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.2"

[[deps.BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

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
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

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

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

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

[[deps.Convex]]
deps = ["AbstractTrees", "BenchmarkTools", "LDLFactorizations", "LinearAlgebra", "MathOptInterface", "OrderedCollections", "SparseArrays", "Test"]
git-tree-sha1 = "145c5e0b3ea3c9dd3bba134a58bab4112aa250c8"
uuid = "f65535da-76fb-5f13-bab9-19810c17039a"
version = "0.14.18"

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
git-tree-sha1 = "cfdfef912b7f93e4b848e80b9befdf9e331bc05a"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.1"

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

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

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
git-tree-sha1 = "6a8dc9f82e5ce28279b6e3e2cea9421154f5bd0d"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.37"

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

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

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
git-tree-sha1 = "f97acd98255568c3c9b416c5a3cf246c1315771b"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.0+0"

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
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "8d70835a3759cdd75881426fced1508bb7b7e1b6"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.1"

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
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

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

[[deps.LDLFactorizations]]
deps = ["AMD", "LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "399bbe845e06e1c2d44ebb241f554d45eaf66788"
uuid = "40e66cde-538c-5869-a4ad-c39174c6795b"
version = "0.8.1"

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

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

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

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "92b7de61ecb616562fd2501334f729cc9db2a9a6"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.10.6"

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

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "7bb6853d9afec54019c1397c6eb610b9b9a19525"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.3.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

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

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "916077e0f0f8966eb0dc98a5c39921fdb8f49eb4"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.0"

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

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

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
git-tree-sha1 = "68604313ed59f0408313228ba09e79252e4b2da8"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "7eda8e2a61e35b7f553172ef3d9eaa5e4e76d92e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "fed057115644d04fba7f4d768faeeeff6ad11a60"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.27"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

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
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

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
git-tree-sha1 = "ee885e0f773804f046fd43d0d4ace305b3d540e2"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.13"

[[deps.SCS]]
deps = ["BinaryProvider", "Libdl", "LinearAlgebra", "MathOptInterface", "Requires", "SCS_GPU_jll", "SCS_jll", "SparseArrays"]
git-tree-sha1 = "c819d023621358f3c08f08d41bd9354cf1357d35"
uuid = "c946c3f1-0d1f-5ce8-9dea-7daa1f7e2d13"
version = "0.8.1"

[[deps.SCS_GPU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "a96402e3b494a8bbec61b1adb86d4be04112c646"
uuid = "af6e375f-46ec-5fa0-b791-491b0dfa44a4"
version = "2.1.4+0"

[[deps.SCS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "6cdaccb5e6a69455f960de1ae445ba1de5db9d0d"
uuid = "f4f2fc5b-1d94-523c-97ea-2ab488bedf4b"
version = "2.1.2+1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "244586bc07462d22aed0113af9c731f2a518c93e"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.10"

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
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

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
git-tree-sha1 = "e1e5ed9669d5521d4bbdd4fab9f0945a0ffceba2"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.14.30"

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
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

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

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

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
# ╟─89b5c4d0-68cb-499f-bf99-216e38b40ca0
# ╟─489464b4-3139-466c-80c7-84449dcec698
# ╟─b4f60458-0299-4246-b2ad-dc3f1db6382b
# ╟─34bf07e8-5c47-4aa8-aa2c-161709c158be
# ╟─65a0683e-3124-4844-ad0a-cccdc9192d05
# ╟─e603c785-b697-48e3-92be-e9213e72215b
# ╟─badff812-6b68-4f6d-b495-3f344873d45c
# ╟─6a4164c2-bb11-437a-bc7d-a12e363d3e84
# ╟─c81486cb-f78e-4eec-a1cd-1ea113428bd6
# ╟─880b1174-0925-4e8c-b9f6-f6e295412824
# ╟─4c2bbb6b-f287-43c4-9e24-f44fb4c87e36
# ╟─00120ff5-e9ed-4e71-8208-cf8efd2eac6a
# ╟─e3ba3cf8-2f94-4276-98ed-6fedcfaadd43
# ╟─a64f7082-b834-4400-880a-032cb9aafe4c
# ╟─a1b50b95-bf7d-4f2a-817d-9edb3418aeb3
# ╟─1852d4c5-e73d-4038-a565-b8fb6ff63502
# ╟─1c9f096f-7dce-4d0c-8a71-6fbe682e514d
# ╟─eb4a5874-c13c-4915-96b9-0c47eaa7f50c
# ╟─66bfd748-3b40-42df-86fa-e55019dba856
# ╟─cf11b553-13e5-488d-bf1c-16de5911d658
# ╟─2353a285-71de-43b2-a60f-5a3274ff9e6b
# ╠═fe2848df-823a-4ed0-918c-2c200957ee80
# ╟─43d75d79-b710-4dc5-9478-dd2b08616be9
# ╟─ebc5ed32-1fe3-4854-a326-7a068e14164b
# ╟─1b25e7d1-909f-4fdc-9700-5d131251e1b5
# ╟─9943d000-83d0-413d-a231-0295fb19df71
# ╟─f66a480b-3f0c-4ebf-a8b8-e0f91dff851d
# ╟─f94ecbae-e823-46f6-a69f-eb7edce7cfe5
# ╟─dc5f6727-cd17-4c75-80ce-94915a0e359a
# ╟─9b9d7fdf-91b1-46e4-bd66-43e50add56be
# ╟─0e7d1741-88a1-4e8e-b964-e8ead4d1807e
# ╟─300b62c8-830a-472f-a07c-17153468c1fb
# ╠═54efa70c-bac6-4d7c-93df-0dfd1b89769d
# ╠═866cc84d-86c1-40e2-bd29-deae01da9a2e
# ╠═7bcbdc4e-a38a-4201-a1ec-d2e4df4d2f6a
# ╟─d1edad45-5df3-43cd-8abc-97d84fab699b
# ╟─bb7582b4-3ecc-44a7-810a-aae0dc2fe816
# ╟─34db9c1e-6af7-4710-8de4-fd0caacbde36
# ╟─46fef026-6bac-401b-ab6f-75f2aff79e6a
# ╟─455b2bea-4d79-4dd8-964c-80835bb88727
# ╟─4d4230e5-28dc-4f87-b732-686a5c91fc4f
# ╠═5c5d5eeb-6775-452f-880d-7b4fa2acda57
# ╟─e2eb17d0-2ca3-437c-97f7-04e76ee879cb
# ╟─a3a08666-f494-489d-a888-8fe7baa70729
# ╟─b05cac4c-d5d9-47f7-ab38-d9be8322ab45
# ╟─61ab2949-d72f-4d80-a717-4b6a9227de0e
# ╠═34b06415-21c1-4904-97f0-ab614447355c
# ╟─e461b560-435b-4104-b42e-af04b1d25984
# ╟─b7b618c2-d48b-45c7-aff8-600eda010169
# ╟─a39b90ec-a4c5-472f-b00e-c81ee9c5576f
# ╟─cbbd8670-49ab-4601-b8d7-9f3f456752e8
# ╟─3f9879fa-ec04-4d13-8207-2c77d9ddfca2
# ╟─22992de2-dfba-4a18-8909-62bddbf4e0e7
# ╟─367158be-9a7f-4b76-96c9-2806f9aad75e
# ╟─fb83174c-fcde-4496-8232-545f17ac9d2d
# ╟─edfbf364-e126-4e95-93d2-a6adfb340045
# ╠═a786ca10-06d2-4b76-97a9-2bcf879ea6cb
# ╠═012ef00f-6176-4d42-801f-5765c47df7ba
# ╟─f067b633-fde0-4743-bd76-f8c390a90950
# ╟─a3d29aa3-96ca-4681-960c-3b4b04b1e40d
# ╟─6bf06c12-cf25-43c4-81f3-b1d79d13fc94
# ╟─4dff7597-bc0e-4288-9413-c26a132c1e44
# ╟─1d72b291-24b7-4ec6-8307-1da0bc4a9183
# ╟─1867ecec-3c3d-4b2b-9036-1488e2184c40
# ╟─0e09d312-2ddf-4d1f-8ad5-a50fb48ca4dd
# ╠═a3f1710e-eb98-46c6-aadb-5cf0b98e1bc6
# ╠═9380908a-8cc5-4d3a-9d6c-03b491198bc1
# ╟─a07d661a-a0b4-4000-b7a4-5f17cda5edc3
# ╠═f0ee3633-1d28-4344-9b85-f5433679582a
# ╠═345d8feb-cfd0-4acb-a626-dc16181ddb68
# ╟─fdc34643-2c27-4bc9-a077-8bac68dc56b7
# ╟─8897de5b-a6c7-4b05-98c6-d738bbb527af
# ╟─379373b1-e563-4341-9978-5b35c768b5c7
# ╟─8b2725a2-8007-46b8-a160-75042562794d
# ╟─95495255-385a-424b-9b84-dcd8a1187282
# ╟─d8af95e4-fe4c-4cf0-8195-dab80629177c
# ╟─3f62e075-4724-456e-ab59-9b85d4263ee6
# ╟─bc4e93f5-7e76-4115-8cfd-039212181fb3
# ╠═20bba789-504b-4d27-a2fc-4f08badc9a53
# ╟─587c70c4-a3b8-4300-b7a4-aa9f6b1bc276
# ╟─421e213e-9780-4a4d-a411-009541c44a9e
# ╟─2d40603b-56e3-49ca-a8ea-e13c933f5e19
# ╟─1c816d04-3b76-4c0c-9ed9-9b20b5ad50d9
# ╟─c39168a1-dee2-4421-8f44-266df47dd08e
# ╟─4e6e394d-0d95-46cc-8344-0222b0fb761e
# ╟─0454a796-b35d-4caa-b847-f578191f2896
# ╟─10fa507e-1429-4eb0-b74c-1e6638725690
# ╟─92673e3f-e436-4c7e-accf-216574233d58
# ╟─e2df85a5-975c-472b-8538-d07b682d1647
# ╟─2e416e6d-ea7c-464c-9632-83b0f7f06b2a
# ╠═5a3500c2-4f82-43e9-a31b-d530f56fdbe9
# ╟─3d1a1f06-5707-46ff-b44a-539c62fe008b
# ╟─90a0a645-cd70-4edb-8240-a152d6e7bb3a
# ╟─cdd66194-cb79-433e-a16c-64be76de83a4
# ╟─4f8a3476-93a8-414e-b169-7046f1a57547
# ╠═ff11bdd2-8ab6-40c0-844b-a3474d7e9a04
# ╟─21791e69-db82-4081-97ec-e7f2c0a46b5d
# ╟─1f9459dd-57ba-4353-81b2-9e04aa8257af
# ╟─90722894-dc53-4964-8141-e35adfeb8a76
# ╟─4d97d757-ed9f-4e77-9b54-6b9c8f2f49ee
# ╟─9a3d5138-0490-42db-9f99-310d94124951
# ╟─a1e1d5f8-e06e-4682-ab54-a9454a8e3b30
# ╟─e36979d5-c1b6-4c17-a65a-d8de8e6bd8d0
# ╟─aeafe1ed-f217-48fd-9624-add5f6f791e6
# ╟─bdf1e4d0-d2e6-4ab5-92c5-5f25518a9acc
# ╟─cac02388-31c4-40b9-8288-a6685e1854fb
# ╟─0a55c3a1-a834-4ec4-beac-0788091f7a70
# ╟─4725bf3c-f058-4238-8c8f-297dc2851c6e
# ╟─b547311c-ddf0-4053-9de4-f0e85b861e63
# ╟─2ed2f3c3-619b-4aed-b88b-b92a43578d84
# ╟─3b9d8429-375d-46c3-a115-eb9022262e09
# ╟─23ae5b4c-f690-46a6-b2c9-3d753cf247d6
# ╟─edf426b6-a571-4938-890a-01a089d02b29
# ╟─c32725a4-e276-4372-8d06-d40ba52c9f09
# ╟─b04cba56-dd48-403b-82fc-1cf3713853a7
# ╟─f4bd98b5-f4a8-424a-b974-41aceede92fb
# ╟─6086ed53-4fbe-4037-b435-8d8aa20a1417
# ╟─defa9be8-4d19-488b-820c-2a3526401cbf
# ╠═f1a71f47-fb19-4988-a439-2ff8d38be5b7
# ╠═91dae79f-e454-4b27-84a7-4cbc6bc33265
# ╠═c8c3fe32-560d-11ec-0617-2dc33608384a
# ╠═5394b13a-e629-47c8-902d-685c061b37ae
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
