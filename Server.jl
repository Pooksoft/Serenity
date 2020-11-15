# include -
include("Include.jl")

# setup shared resources -
const SERENITY_ROUTER = HTTP.Router()

# Routes -
function echo(req::HTTP.Request)
    return HTTP.Response(200, "Echo")
end

# Register routes, and start the server -
HTTP.@register(SERENITY_ROUTER,"GET","/pooksoft/serenity/api/v1/echo", echo)
HTTP.serve(SERENITY_ROUTER, Sockets.localhost, 8000)