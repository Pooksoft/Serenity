# fire up the server -
using HTTP
using JSON

# create a dict -
json_body = Dict{String,Any}()
json_body["user_email_address"] = "jvarner@pooksoft.com"
json_body["user_api_key"] = "a404b4f0-e54a-44cb-9759-636338efc3c4"

# go -
raw_response = HTTP.request("POST","http://localhost:8000/pooksoft/serenity/api/v1/authenticate",
    ["Content-Type"=>"application/json"], 
    JSON.json(json_body))

dd = JSON.parse(String(raw_response.body))
