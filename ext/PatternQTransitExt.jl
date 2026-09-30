# Optional transit response formats (query(...; format="transit+json" /
# "transit+msgpack")), loaded when Transit.jl is installed and imported.
# Decoded values are converted to what the JSON path gives (tonative of
# JSON3): OrderedDict{String,Any} maps, Vector{Any} arrays, ":ns/name"
# strings for keywords, the service's ISO strings for instants.
module PatternQTransitExt

using Dates
using OrderedCollections
import PatternQ
import Transit

plain(x::Union{String,Bool,Nothing,Int64,Float64}) = x
plain(x::AbstractFloat) = Float64(x)          # includes Decimal
plain(x::Integer) = x
plain(x::Symbol) = ":" * String(x)             # keywords
plain(x::Transit.TSymbol) = ":" * x.s
plain(x::Base.UUID) = string(x)
plain(x::Transit.TaggedValue) = plain(x.value)
plain(x::AbstractString) = String(x)
function plain(x::DateTime)
    # as the JSON response writes instants: milliseconds without trailing zeros
    s = Dates.format(x, dateformat"yyyy-mm-ddTHH:MM:SS")
    frac = rstrip(lpad(string(Dates.millisecond(x)), 3, '0'), '0')
    isempty(frac) ? s * "Z" : s * "." * frac * "Z"
end
function plain(x::AbstractVector)
    out = Vector{Any}(undef, length(x))
    @inbounds for i in eachindex(x)
        v = x[i]
        out[i] = v isa Union{String,Int64,Float64,Bool,Nothing} ? v : plain(v)
    end
    out
end
plain(x::AbstractDict) = OrderedDict{String,Any}(plain_key(k) => plain(v) for (k, v) in x)
plain(x) = Any[plain(v) for v in x]            # sets, lists and other collections

plain_key(k::String) = k
plain_key(k) = string(plain(k))

function PatternQ.decode_transit(body::Vector{UInt8}, format::AbstractString)
    fmt = format == "transit+msgpack" ? :msgpack : :json
    plain(Transit.parse(IOBuffer(body); format=fmt))
end

end
