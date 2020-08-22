# exports

export GeneralBC, BCTypes, PeriodicBoundary, ReflectiveBoundary, OpenBoundary
export CubicPBC, apply_simulation_bc!, applybc

abstract type SimulationBoundaries end
abstract type BCTypes end

struct PeriodicBoundary{T <: Real} <: BCTypes
    X0::T
    X1::T
    L::T
end

function PeriodicBoundary(X0::Real, X1::Real)
    L = X1 - X0
    return PeriodicBoundary(promote(X0, X1, L)...)
end

struct ReflectiveBoundary{T <: Real} <: BCTypes
    X0::T
    X1::T
    L::T
    skin::T
end

function ReflectiveBoundary(X0::Real, X1::Real; skin=0.001)
    L = X1 - X0
    X0 += skin*L
    X1 -= skin*L
    L = X1 - X0
    return ReflectiveBoundary(promote(X0, X1, L, skin)...)
end

struct OpenBoundary <: BCTypes
end

struct GeneralBC{T1 <: BCTypes, T2 <: BCTypes, T3 <: BCTypes} <: SimulationBoundaries
    X::T1
    Y::T2
    Z::T3
end

function CubicPBC(L::Real)
    X = PeriodicBoundary(0.0,L,L)
    GeneralBC(X,X,X)
end

function apply_simulation_bc!(x, v, BC::GeneralBC)
    for i in 1:size(x)[2]
        x[1,i], v[1,i] = applybc(x[1,i], v[1,i], BC.X)
        x[2,i], v[2,i] = applybc(x[2,i], v[2,i], BC.X)
        x[3,i], v[3,i] = applybc(x[3,i], v[3,i], BC.X)
    end
end

function applybc(x, v, BC::OpenBoundary)
    return x, v
end

function applybc(x, v, BC::PeriodicBoundary)
    return x - BC.L*floor((x-BC.X0)/BC.L), v
end

function applybc(x, v, BC::ReflectiveBoundary)
    if x>BC.X1 || x<BC.X0
        return x, -v
    else
        return x, v
    end
end

function delta(dX, BC::Union{ReflectiveBoundary,OpenBoundary})
    return dX
end

function delta(dX, BC::PeriodicBoundary)
    abs_dX = abs(dX)
    if 2*abs_dX > BC.L
        abs_dX = BC.L-abs_dX
        dX = -abs_dX*dX/abs(dX)
    end
    dX
end


function distance(a, b, BC)
    dX = delta(a[1]-b[1],BC.X)
    dY = delta(a[2]-b[2],BC.Y)
    dZ = delta(a[3]-b[3],BC.Z)
    r2 = dX*dX+dY*dY+dZ*dZ
    r = sqrt(r2)
    dr = [dX, dY, dZ]
    return r, r2, dr
end



#
