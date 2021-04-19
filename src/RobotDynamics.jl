module RobotDynamics

using Rotations
using StaticArrays
using LinearAlgebra
using ForwardDiff
using UnsafeArrays
using RecipesBase

using Rotations: skew
using StaticArrays: SUnitRange

# include("rbstate.jl")
# include("jacobian.jl")
# include("knotpoint.jl")
include("model.jl")
# include("liestate.jl")
# include("rigidbody.jl")
# include("integration.jl")
# include("trajectories.jl")
# include("linearmodel.jl")
# include("linearization.jl")
# include("plot_recipes.jl")

export
    AbstractModel,
    discrete_dynamics,
    discrete_dynamics!,
    discrete_jacobian!,
    # DynamicsExpansion,
    # dynamics,
    # jacobian!,

    # linearize,
    # linearize!,
    # state_dim,
    # control_dim,
    # state_diff_size,
    # rollout!

# # rigid bodies
# export
#     LieGroupModel,
#     RigidBody,
#     RBState,
#     orientation,
#     linear_velocity,
#     angular_velocity

# # linear model
# export
#     LinearModel,
#     linear_dynamics,
#     LinearizedModel,
#     linearize_and_discretize!,
#     discretize,
#     discretize!,
#     update_trajectory!

# knotpoints
# export
#     AbstractKnotPoint,
#     KnotPoint,
#     StaticKnotPoint,
#     Traj,
#     state,
#     control,
#     states,
#     controls,
#     set_states!,
#     set_controls!

# integration
export
    QuadratureRule
    # RK2,
    # RK3,
    # RK4,
    # HermiteSimpson,
    # PassThrough,
    # Exponential

end # module
