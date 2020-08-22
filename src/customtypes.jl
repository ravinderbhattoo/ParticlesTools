# exports
export WellArray, fillit!, resetit!

mutable struct WellArray{T <: Real}
    val::Array{T,1}
    ind::Int64
end

function WellArray(N; fill_with=0.0)
    c = Array{Float64}(undef, N)
    fill!(c, fill_with)
    WellArray(c, 0)
end

function fillit!(m::WellArray, val)
    m.ind += 1
    m.val[m.ind] = val
    nothing
end

function resetit!(m::WellArray; fill_with=0.0)
    m.ind = 0
    fill!(m.val, fill_with)
end


Base.iterate(m::WellArray, state=1) = state > m.ind  ? nothing : (m.val[state], state+1)

Base.length(m::WellArray) = m.ind
Base.lastindex(m::WellArray) = m.ind
Base.firstindex(m::WellArray) = 1

Base.getindex(m::WellArray, i::Int64) = i > m.ind ? throw(BoundsError([m],["$i with length $(m.ind)"])) : getindex(m.val, i)

Base.getindex(m::WellArray{S}, I::AbstractRange{Int}) where S =  last(I) > m.ind ? throw(BoundsError([m],["$I with length $(m.ind)"])) : S[ m.val[i] for i in I ]

Base.show(stream::IO, m::WellArray{S}) where S = Base.show(stream, m.val[1:m.ind])

Base.:+(m1::WellArray{S}, m2::WellArray{S}) where S = Base.:+(m1[1:end], m2[1:end])

Base.:-(m1::WellArray{S}, m2::WellArray{S}) where S = Base.:-(m1[1:end], m2[1:end])

#
