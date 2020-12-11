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

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == false)

        # build the payload dictionary -
        payload_dictionary = JSON.parse(body_string)

        # grab the password -
        username = payload_dictionary["payload"]["username"]
        password = payload_dictionary["payload"]["password"]
        if (password == "1234")

            # create a UUID -
            uuid_username = string(UUIDs.uuid4())

            # setup the response -
            response_json_dictionary["usertoken"] = uuid_username
            response_json_dictionary["authenticated_status"] = 1

            # cache the uuid -
            SERENITY_SESSION[username] = uuid_username
        else
            response_json_dictionary["authenticated_status"] = 0
        end
    else
        response_json_dictionary["authenticated_status"] = 0
    end

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

function expiration(request::HTTP.Request)

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # setup -
    response_json_dictionary["authenticated_status"] = 0

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == false)

        # build the payload dictionary -
        payload_dictionary = JSON.parse(body_string)

        # grab the username, and token - 
        # we need to check  this against the cached usertoken -
        username = payload_dictionary["user"]["username"]
        usertoken = payload_dictionary["user"]["usertoken"]
        
        # ok, do we have this username key?
        if (haskey(SERENITY_SESSION,username) == true)
            
            # check: r we cached?
            cached_user_token = SERENITY_SESSION[username]
            if (usertoken == cached_user_token)
                response_json_dictionary["authenticated_status"] = 1

                # ok: so we good to go re authentication -
                # ...

            end
        end
    end

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

# Register routes, and start the server -
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/expiration", expiration)
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/authenticate", authenticate)
HTTP.@register(SERENITY_ROUTER,"GET","/pooksoft/serenity/api/v1/echo", echo)
HTTP.serve(SERENITY_ROUTER, Sockets.localhost, 8000)