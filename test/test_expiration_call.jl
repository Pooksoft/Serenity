# fire up the server -
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

# -- STEP 2: COMPUTE CALL ---------------------------------------------- #
# initialize -
json_body = Dict{String,Any}()

# setup user authentication -
json_body["user_api_key"] = "a404b4f0-e54a-44cb-9759-636338efc3c4"
json_body["user_session_token"] = user_session_token

# load the config -
path_to_config_file = "./test/config/Straddle.json"
contract_dictionary = JSON.parsefile(path_to_config_file)
json_body["contract_set"] = contract_dictionary

# set the price range -
simulation_dictionary = Dict{String,Any}()
simulation_dictionary["start_price"] = 135.0
simulation_dictionary["stop_price"] = 174.0
simulation_dictionary["length"] = 100
json_body["simulated_price_range"] = simulation_dictionary

# go -
raw_response = HTTP.request("POST","http://localhost:8000/pooksoft/serenity/api/v1/contract/expiration",
    ["Content-Type"=>"application/json"], 
    JSON.json(json_body))

# get dictionary back from call -
dd = JSON.parse(String(raw_response.body))

# pull out the compute_result_array -
compute_result_array = dd["compute_result_array"]
price_array = compute_result_array[1]
profit_array = compute_result_array[2]
plot(price_array, profit_array, lw=2)
xlabel!("Share Price (USD/share)",fontsize=14)
ylabel!("Profit/Loss (USD/share)",fontsize=14)
# ---------------------------------------------------------------------- #