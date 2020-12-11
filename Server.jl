# include -
include("Include.jl")

# setup shared resources -
const SERENITY_ROUTER = HTTP.Router()
const SERENITY_SESSION = Dict{String,Any}()

# setup paths -
const path_to_server = "/Users/jeffreyvarner/Desktop/julia_work/Serenity" 
const path_to_database_file = "$(path_to_server)/database/Serenity.db"
const SERENITY_DB_CONNECTION = SQLite.DB(path_to_database_file)

# Routes -
"""
    echo()
"""
function echo(req::HTTP.Request)
    return HTTP.Response(200, "Echo")
end

"""
    authenticate()
"""
function authenticate(request::HTTP.Request)

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # default is *NOT* authenticated -
    response_json_dictionary["user_authenticated_status"] = 0

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == false)

        # build the payload dictionary -
        payload_dictionary = JSON.parse(body_string)

        # grab the user information -
        username = payload_dictionary["user_email_address"]
        password = payload_dictionary["user_api_key"]

        # build a PSAuthenticationUserModel
        user_authentication_model = PSAuthenticationUserModel(username,password)
        authentication_result = authenticate_user_api_call(SERENITY_DB_CONNECTION, user_authentication_model)
        if (isa(authentication_result.value, Exception) == true)
            # what?
        end
        authentication_status = authentication_result.value

        # status = true?
        if (authentication_status == true)

            # create a UUID -
            user_session_token = string(UUIDs.uuid4())

            # setup the response -
            response_json_dictionary["user_session_token"] = user_session_token
            response_json_dictionary["user_authenticated_status"] = 1

            # cache the uuid -
            SERENITY_SESSION[username] = user_session_token
        end        
    end

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

function compute_contract_set_expiration(request::HTTP.Request)

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == false)

        # build the payload dictionary -
        payload_dictionary = JSON.parse(body_string)

        # call -
        compute_result = compute_contract_set_expiration(payload_dictionary)
        if (isa(compute_result.value, Exception) == true)
            
            # what error?
            # ...

        end
        pla = compute_result.value

        # get the size of the output -
        (number_of_timesteps,number_of_cols) = size(pla)
        for timestep_index = 1:number_of_timesteps
        end
    end

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

# Register routes, and start the server -
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/expiration", compute_contract_set_expiration)
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/authenticate", authenticate)
HTTP.@register(SERENITY_ROUTER,"GET","/pooksoft/serenity/api/v1/echo", echo)
HTTP.serve(SERENITY_ROUTER, Sockets.localhost, 8000)