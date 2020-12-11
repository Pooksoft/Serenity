struct PSAuthenticationUserModel

    # data -
    user_email_address::String
    user_api_key::String

    function PSAuthenticationUserModel(user_email_address::String, user_api_key::String)
        this = new(user_email_address,user_api_key)
    end
end