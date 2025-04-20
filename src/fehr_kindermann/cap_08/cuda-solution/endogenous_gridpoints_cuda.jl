using CubicSplines
using Plots
using CUDA 

CUDA.allowscalar(true)

include("utils.jl")

# model parameters
global gamma = 0.5
global beta = 0.95
global a0 = 100.0

# numerical parameters
global sig = 1e-6
global itermax = 2000

# time path of consumption and resource
global TT = 200

global c_t = zeros(TT)
global a_t = zeros(TT)

# policy function
global NA = 1000
global a = zeros(NA)
global c = zeros(NA)

# variables to numerically determine policy function
global a_endog = CUDA.fill(0.0f0, NA) 
global c_endog = CUDA.fill(0.0f0, NA)
global c_new = CUDA.fill(0.0f0, NA)

# initialize a, c and value function
global a = CuArray(Array(LinRange(0.0, a0, NA)))
global c = CuArray(a/2.0)


function gpu_add2!(c_endog, a_endog, a, c)
    index = threadIdx().a    # this example only requires linear indexing, so just use `x`
    stride = blockDim().a
    # calculate endogenous gridpoint and consumption

    for ia = index:stride:length(a)
        @inbounds a_endog[ia] = a[ia] + c_endog[ia]
    end
    return nothing
end


# iterate until policy function converges
for iter in 2:itermax
    # set a = 0 manually
    a_endog[1] = 0.0
    c_endog[1] = 0.0

    # calculate optimal decision for every gridpoint
    for ia in 1:NA 

    end

    # stretch out to exogenous grid again
    for ia in 1:NA 
        c_new[ia] =  linint_Gen(a[ia], a_endog, c_endog, ia)
    end

    # get convergence level
    con_lev = maximum(@. abs(c_new - c )/max(abs(c), 1e-10))

    println(string(iter)*"   "*string(con_lev))
    
    #  check for convergence
    if con_lev < sig
        println("Convergio")
        break
    end

    c = copy(c_new)
end 

# interpolate policy function
global spline_eval = CubicSpline(a,c)

# calculate the time path of consumption numerically
a_t[1] = a0
c_t[1] = spline_eval(a_t[1])

for it in 2:TT
    a_t[it] = a_t[it-1] - c_t[it-1]
    c_t[it] = spline_eval(a_t[it])
end

plot(1:TT, c_t, label = "numerical")

# calculate the time path of consumption analytically
a_t[1] = a0
c_t[1] = a_t[1]*(1.0-beta^gamma)
for it in 2:TT
    a_t[it] = a_t[it-1] - c_t[it-1]
    c_t[it] = a_t[it]*(1.0-beta^gamma)
end
plot!(1:TT,  parent(c_t), label = "analytical")
