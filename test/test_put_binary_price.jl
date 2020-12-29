# we need these packages ...
using HTTP
using JSON
using Plots

# -- STEP 1: AUTHENTICATION --------------------------------------------- #
# need to authenticate first -
# create a dict -
json_authentication_body = Dict{String,Any}()
json_authentication_body["user_email_address"] = "jvarner@pooksoft.com"
json_authentication_body["user_api_key"] = "a404b4f0-e54a-44cb-9759-636338efc3c4"

# go for authentication ...
raw_response = HTTP.request("POST","http://localhost:8000/pooksoft/serenity/api/v1/authenticate",
    ["Content-Type"=>"application/json"], 
    JSON.json(json_authentication_body))

# if this worked, this should have the session token
authentication_dictionary = JSON.parse(String(raw_response.body))
user_session_token = authentication_dictionary["user_session_token"]
# ---------------------------------------------------------------------- #

# -- STEP 2: COMPUTE --------------------------------------------------- #
# initialize -
json_body = Dict{String,Any}()

# setup user authentication -
json_body["user_api_key"] = "a404b4f0-e54a-44cb-9759-636338efc3c4"
json_body["user_session_token"] = user_session_token

# load the example data -
path_to_contract_file = "./test/config/Put-Binary-ALLY.json"
contract_dictionary = JSON.parsefile(path_to_contract_file)
json_body["contract_set_parameters"] = contract_dictionary["contract_set_parameters"]
json_body["binary_lattice_parameters"] = contract_dictionary["binary_simulation_parameters"]

# set the underlying price?
baseUnderlyingPrice = 34.54
json_body["underlying_price_value"] = baseUnderlyingPrice

# setup the strike price range -
simulation_dictionary = Dict{String,Any}()
simulation_dictionary["start_price"] = 28.0
simulation_dictionary["stop_price"] = 45.0
simulation_dictionary["number_of_steps"] = 200
json_body["simulated_strike_price_dictionary"] = simulation_dictionary

# go -
raw_response = HTTP.request("POST","http://localhost:8000/pooksoft/serenity/api/v1/contract/binary/put/price",
    ["Content-Type"=>"application/json"], 
    JSON.json(json_body))

# get dictionary back from call -
dd = JSON.parse(String(raw_response.body))
# ---------------------------------------------------------------------- #

# -- PLOT -------------------------------------------------------------- #
pa = dd["compute_result_array"]
spa = pa[1]
iva = pa[3]
opa = pa[4]
# ---------------------------------------------------------------------- #