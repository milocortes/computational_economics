using CubicSplines
using Plots
using Roots

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
global c_new = zeros(NA)

# variables to communicate with function
global ia_com

# the first order condition
function foc(x_in, a_val)
    # calculate right hand side of foc
    cplus = spline_eval(x_in)

    # get foc
    return a_val - x_in - beta^(-gamma)*cplus
end

# initialize a and c
global a = Array(LinRange(0.0, a0, NA))
global c = a/2.0

global spline_eval = CubicSpline(a,c)

# iterate until policy function converges
@elapsed for iter in 2:itermax
    # set a = 0 manually
    c_new[1] = 0.0

    # calculate optimal decision for every gridpoint
    Threads.@threads for ia in 1:NA
        # initialize starting value and communicate resources
        x_in = a[ia] - c[ia]        
        # solve the household problem using rootfinding
        x_root = find_zero(foc, x_in, xatol=1e-12, p = a[ia])
        # get optimal consumption and value function
        c_new[ia] = a[ia] - x_root
    end

    # interpolate coefficients
    spline_eval = CubicSpline(a,c_new)

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
