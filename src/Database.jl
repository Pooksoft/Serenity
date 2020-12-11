function authenticate_user_api_call(db::SQLite.DB, model::PSAuthenticationUserModel)::PooksoftBase.PSResult

    # get data from the user model -
    user_email_address = model.user_email_address
    test_user_api_key = model.user_api_key

    # formulate the SQL call -
    sql_string = "SELECT user_api_key FROM USER_REGISTRATION_TABLE WHERE user_email_address='$(user_email_address)';"

    # execute the call -
    query_result = DBInterface.execute(db,sql_string)

    # turn the result into a table (and return for now so we can see what gets returned)
    df = DataFrame(query_result)
    
    # case 1: no records = ok, so we don't have any records, return false
    if (isempty(df) == true)
        return PSResult{Bool}(false)
    end
    
    # case 2: we have records, but more than one?
    if (size(df,1)>1)
        return PSResult{Bool}(false) 
    end

    # case 3: we have 1 record, check: do the key's match?
    stored_user_api_key = df[!,Symbol("user_api_key")][1]
    if (stored_user_api_key != test_user_api_key)
        return PSResult{Bool}(false)
    end

    # return -
    return PSResult{Bool}(true)
end