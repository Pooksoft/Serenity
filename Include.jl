# Setup the paths -
const _PATH_TO_ROOT = pwd()
const _PATH_TO_SRC = "$(_PATH_TO_ROOT)/src"
const _PATH_TO_DATA = "$(_PATH_TO_ROOT)/data"
const _PATH_TO_CONFIG = "$(_PATH_TO_ROOT)/config"

# activate the project -
import Pkg
Pkg.activate(_PATH_TO_ROOT);

# load external packages -
using HTTP
using JSON
using Sockets
using UUIDs
using Logging
using SQLite
using DataFrames
using Query
using CSV
using Plots
using Dates

# Pooksoft codes -
using PooksoftBase
using PooksoftAlphaVantageDataStore
using PooksoftAssetModelingKit
using PooksoftOptionsKit

# include my codes -
include("$(_PATH_TO_SRC)/Types.jl")
include("$(_PATH_TO_SRC)/Database.jl")
include("$(_PATH_TO_SRC)/Options.jl")
include("$(_PATH_TO_SRC)/Model.jl")
include("$(_PATH_TO_SRC)/Analysis.jl")
include("$(_PATH_TO_SRC)/Checks.jl")

# notify the user
@info "Starting the Serenity server: http://localhost:8000"