# fire up the server -
using HTTP
using JSON

# initialize -
json_body = Dict{String,Any}()

# setup user authentication -
json_body["user_email_address"] = "jvarner@pooksoft.com"
json_body["user_api_key"] = "a404b4f0-e54a-44cb-9759-636338efc3c4"

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
raw_response = HTTP.request("POST","http://localhost:8000/pooksoft/serenity/api/v1/expiration",
    ["Content-Type"=>"application/json"], 
    JSON.json(json_body))

dd = JSON.parse(String(raw_response.body))

