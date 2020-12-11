include("../Include.jl")

# setup paths -
const path_to_server = "/Users/jeffreyvarner/Desktop/julia_work/Serenity" 
const path_to_database_file = "$(path_to_server)/database/Serenity.db"
const SERENITY_DB_CONNECTION = SQLite.DB(path_to_database_file)

# Generate a user model -
user_authentication_model = PSAuthenticationUserModel("jvarner@pooksoft.com","a404b4f0-e54a-44cb-9759-636338efc3c4")

# execute the db call -
result = authenticate_user_api_call(SERENITY_DB_CONNECTION, user_authentication_model)
if (isa(result.value,Exception) == true)
    throw(result.value)
end
df = result.value
