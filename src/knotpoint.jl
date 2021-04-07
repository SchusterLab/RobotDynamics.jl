export
    GeneralKnotPoint,
    state,
    control


"""
	AbstractKnotPoint

Stores the states, controls, time, and time step at a single point along a trajectory of a
forced dynamical system.

# Interface
All instances of `AbstractKnotPoint` should support the following methods:

	state(z)     # state vector
	control(z)   # control vector
	z.t::Real    # time
	z.dt::Real   # time to next point (time step)

Alternatively, the methods `state` and `control` will be automatically defined if the
following fields are present:
- `z.z`: the stacked vector `[x;u]`
- `z.ix`: the indices of the states, such that `x = z.z[z.ix]`
- `z.iu`: the indices of the controls, such that `x = u.z[z.iu]`
"""
abstract type AbstractKnotPoint end


"""
	state(::AbstractKnotPoint)

Return the `n`-dimensional state vector
"""
@inline state(z::AbstractKnotPoint) = z.z[z.ix]

"""
	control(::AbstractKnotPoint)

Return the `m`-dimensional control vector
"""
@inline control(z::AbstractKnotPoint) = z.z[z.iu]


"""
	is_terminal(::AbstractKnotPoint)::Bool

Determine if the knot point is the terminal knot point.
"""
@inline is_terminal(z::AbstractKnotPoint) = z.terminal


"""
	get_z(::AbstractKnotPoint)

Returns the stacked state-control vector `z`, or just the state vector if `is_terminal(z) == true`.
"""
@inline get_z(z::AbstractKnotPoint) = is_terminal(z) ? state(z) : z.z


"""
	set_state!(z::AbstractKnotPoint, x::AbstractVector)

Set the state in `z` to `x`.
"""
set_state!(z::AbstractKnotPoint, x) = for i in z.ix; z.z[i] = x[i]; end


"""
	set_control!(z::AbstractKnotPoint, u::AbstractVector)

Set the controls in `z` to `u`.
"""
set_control!(z::AbstractKnotPoint, u) = for (i,j) in enumerate(z.iu); z.z[j] = u[i]; end


"""
	set_z!(z::AbstractKnotPoint, z_::AbstractVector)

Set both the states and controls in `z` from the stacked state-control vector `z_`, unless
`is_terminal(z)`, in which case `z_` is assumed to be the terminal states.
"""
set_z!(z::AbstractKnotPoint, z_) = is_terminal(z) ? set_state!(z, z_) : copyto!(z.z, z_)


function Base.isapprox(z1::AbstractKnotPoint, z2::AbstractKnotPoint)
    get_z(z1) ≈ get_z(z2) && z1.t ≈ z2.t && z1.dt ≈ z2.dt
end


function Base.:+(z1::AbstractKnotPoint, z2::AbstractKnotPoint)
	StaticKnotPoint(z1.z + z2.z, z1.ix, z1.iu, z1.dt, z1.t, z1.terminal)
end


function Base.:+(z::Tk, x::AbstractVector) where {Tk<:AbstractKnotPoint}
	Tk(z.z + x, z.ix, z.iu, z.dt, z.t, z1.terminal)
end


@inline Base.:+(x::AbstractVector, z1::AbstractKnotPoint) = z1 + x

function Base.:*(a::Real, z::Tk) where {Tk<:AbstractKnotPoint}
	Tk(z.z*a, z.ix, z.iu, z.dt, z.t, z.terminal)
end


@inline Base.:*(z::AbstractKnotPoint, a::Real) = a*z


"""
	GeneralKnotPoint{Tz,Tix,Tiu,T} <: AbstractKnotPoint

A mutable instantiation of the `AbstractKnotPoint` interface where
the joint vector `z = [x;u]` has type `Tz`, the state indices `ix`
have type `Tix`, the control indices `iu` have type `Tiu`, and `dt` and `t`
have type `T`.
"""
mutable struct GeneralKnotPoint{Tz,Tix,Tiu,T} <: AbstractKnotPoint
    z::Tz
    ix::Tix
    iu::Tiu
    dt::T
    t::T
    terminal::Bool
end

function GeneralKnotPoint(z::Tz, ix::Tix, iu::Tiu, dt::T, t::T, terminal=false) where {Tz,Tix,Tiu,T}
    return GeneralKnotPoint{Tz,Tix,Tiu,T}(z, ix, iu, dt, t, terminal)
end

function Base.copy(z::GeneralKnotPoint)
    GeneralKnotPoint(Base.copy(z.z), z.ix, z.iu, z.dt, z.t, z.terminal)
end

