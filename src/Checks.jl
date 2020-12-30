function check_authentication_status(request_body_dictionary::Dict{String,Any})::Union{Nothing,HTTP.Response}

    # initialize -
    response_json_dictionary = Dict{String,Any}()
    response_code = 200

    # check: do we have a user_api_key?
    if (haskey(request_body_dictionary,"user_api_key") == false)
        
        # ok, so we are missing the user_api_key, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The user_api_key is missing from the request_body_dictionary"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end
    user_api_key = request_body_dictionary["user_api_key"]

    # check: do we have the user session token?
    if (haskey(request_body_dictionary,"user_session_token") == false)
        
        # ok, so we are missing the user_session_token, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The user_session_token is missing from the request_body_dictionary"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end
    user_session_token = request_body_dictionary["user_session_token"]

    # check: does the SERENITY_SESSION have an entry for this user_api_key?
    if (haskey(SERENITY_SESSION, user_api_key) == false)
        
        # ok, so we are missing the user_api_key in the SERENITY_SESSION, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "The user_api_key is missing from SERENITY_SESSION"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end    

    # ok, check the session token, is this a legit token that is matched to this user?
    if (SERENITY_SESSION[user_api_key] != user_session_token)
        
        # ok, so we are missing the user_api_key in the SERENITY_SESSION, send back 
        # the appropriate response code, stack overflow says 422?
        response_code = 422

        # message -
        response_json_dictionary["error_location"] = "$(@__LINE__) of $(@__FILE__)"
        response_json_dictionary["error_message"] = "incorrect user_session_token"

        # encode -
        buffer = JSON.json(response_json_dictionary)
        response = HTTP.Response(response_code, buffer)

        # return -
        return response
    end

    # return -
    return nothing
end