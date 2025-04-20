using CubicSplines
using Plots

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
global NA = 1000000
global a = zeros(NA)
global c = zeros(NA)

# variables to numerically determine policy function
global a_endog = zeros(NA) 
global c_endog = zeros(NA) 
global c_new = zeros(NA) 

# initialize a, c and value function
global a = Array(LinRange(0.0, a0, NA))
global c = a/2.0


# iterate until policy function converges
for iter in 2:itermax
    # set a = 0 manually
    a_endog[1] = 0.0
    c_endog[1] = 0.0

    # calculate optimal decision for every gridpoint
    for ia in 1:NA 
        # calculate endogenous gridpoint and consumption
        c_endog[ia] = c[ia]*beta^(-gamma)
        a_endog[ia] = a[ia] + c_endog[ia]
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
