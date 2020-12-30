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

    # get body of the message -
    body_string = String(request.body)
    if (isempty(body_string) == true)
        
        # ok, so we are missing the body of the message, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so if I get here, then I *have* a body, and its beautiful, so amazing ...
    # build the request_body_dictionary -
    request_body_dictionary = JSON.parse(body_string)    

    # check: do we have the user_api_key and the user_email_address?
    if (haskey(request_body_dictionary,"user_api_key") == false)
        # ok, so we are missing the user_api_key, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    if (haskey(request_body_dictionary,"user_email_address") == false)
        # ok, so we are missing the user_api_key, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so have both the api_key and the user_email address -
    user_email_address = request_body_dictionary["user_email_address"]
    user_api_key = request_body_dictionary["user_api_key"]

    # build a PSAuthenticationUserModel
    user_authentication_model = PSAuthenticationUserModel(user_email_address,user_api_key)
    authentication_result = authenticate_user_api_call(SERENITY_DB_CONNECTION, user_authentication_model)
    if (isa(authentication_result.value, Exception) == true)
        
        # something happend - like the server crashed. 
        # diff than wrong info was passed in the body = we check this below
        # pass back a 501 response code -
        response_code = 501

        # should be provide some error information here?
        # ...

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # check: do we pass the wrong information or some such thing?
    authentication_status = authentication_result.value
    if (authentication_status == false)
        
        # pass back a 101 response code - unathorized
        response_code = 401

        # should be provide some error information here?
        # ...

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so finally - we have an authenticated user.
    # create a UUID -
    user_session_token = string(UUIDs.uuid4())

    # setup the response -
    response_json_dictionary["user_session_token"] = user_session_token
    response_json_dictionary["user_authenticated_status"] = 1

    # cache the uuid -
    SERENITY_SESSION[user_api_key] = user_session_token

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

"""
    compute_put_price_binary_model()
"""
function compute_put_price_binary_model(request::HTTP.Request)

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # get body of the message -
    body_string = String(request.body)
    if (isempty(body_string) == true)
        
        # ok, so we are missing the body of the message, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # error message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The request.body is empty?"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so if I get here, then I *have* a body, and its beautiful, so amazing ...
    # build the request_body_dictionary -
    request_body_dictionary = JSON.parse(body_string)    

    # check the authentication status. This method returns a HTTP.response or nothing
    authenticate_result = check_authentication_status(request_body_dictionary)
    if (isnothing(authenticate_result) == false)
        return authenticate_result
    end
    
    # ok, so I should have everything I need. Do the computation, 
    compute_result = compute_put_price_binary_model(request_body_dictionary)
    if (isa(compute_result.value, Exception) == true)
        
        # ok, so we are missing the user_api_key in the SERENITY_SESSION, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 500

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "Compute failure. Error returned: $(compute_result.value)"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end
    
    # encode -
    response_json_dictionary["result"] = compute_result.value

    # encode the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

"""
    compute_option_contract_set_value_at_expiration()
"""
function compute_option_contract_set_value_at_expiration(request::HTTP.Request)::HTTP.Response

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == true)
        
        # ok, so we are missing the body of the message, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # error message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The request.body is empty?"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so if I get here, then I *have* a body, and its beautiful, so amazing ...
    # build the request_body_dictionary -
    request_body_dictionary = JSON.parse(body_string)    

    # check the authentication status. This method returns a HTTP.response or nothing
    authenticate_result = check_authentication_status(request_body_dictionary)
    if (isnothing(authenticate_result) == false)
        return authenticate_result
    end
    
    # call -
    compute_result = compute_option_contract_set_expiration(request_body_dictionary)
    if (isa(compute_result.value, Exception) == true)
        
        # ok, so we are missing the user_api_key in the SERENITY_SESSION, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 500

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "Compute failure. Error returned: $(compute_result.value)"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end
    pla = compute_result.value
    
    # encode -
    response_json_dictionary["result"] = pla

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

"""
    compute_equity_price_binary_model()
"""
function compute_equity_price_binary_model(request::HTTP.Request)::HTTP.Response

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # check - do we have data from the body?
    body_string = String(request.body)
    if (isempty(body_string) == true)
        
        # ok, so we are missing the body of the message, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # error message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The request.body is empty?"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # ok, so if I get here, then I *have* a body, and its beautiful, so amazing ...
    # build the request_body_dictionary -
    request_body_dictionary = JSON.parse(body_string)    

    # check the authentication status. This method returns a HTTP.response or nothing
    authenticate_result = check_authentication_status(request_body_dictionary)
    if (isnothing(authenticate_result) == false)
        return authenticate_result
    end

    # Compute -
    compute_result = compute_equity_price_binary_model(request_body_dictionary)
    if (isa(compute_result.value, Exception) == true)
        
        # ok, so we are missing the user_api_key in the SERENITY_SESSION, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 500

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "Compute failure. Error returned: $(compute_result.value)"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end
    bpa = compute_result.value

    # encode -
    response_json_dictionary["result"] = bpa

    # grab the payload -
    buffer = JSON.json(response_json_dictionary)
    response = HTTP.Response(response_code, buffer)

    # return -
    return response
end

# Register routes, and start the server -
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/contract/expiration", compute_option_contract_set_value_at_expiration)
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/authenticate", authenticate)
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/contract/binary/put/price", compute_put_price_binary_model)
HTTP.@register(SERENITY_ROUTER, "POST", "/pooksoft/serenity/api/v1/equity/binary/price", compute_equity_price_binary_model)
HTTP.@register(SERENITY_ROUTER,"GET","/pooksoft/serenity/api/v1/echo", echo)
HTTP.serve(SERENITY_ROUTER, Sockets.localhost, 8000)