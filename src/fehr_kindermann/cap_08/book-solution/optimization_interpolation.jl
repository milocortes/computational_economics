using OffsetArrays
using Interpolations
using Optim
using CubicSplines
using Plots


# model parameters
global gamma = 0.5
global egam = 1.0-1.0/gamma
global beta = 0.95
global a0 = 100.0

# numerical parameters
global sig = 1e-6
global itermax = 2000

# time path of consumption and resource
global TT = 200
#global c_t = OffsetArray(zeros(TT+1), 0:TT);
#global a_t = OffsetArray(zeros(TT+1), 0:TT);
global c_t = zeros(TT)
global a_t = zeros(TT)

# value and policy function
global NA = 1000

#global V = OffsetArray(zeros(NA+1), 0:NA)
global V = zeros(NA)

# variables to numerically determine value and policy function
#global V_new = OffsetArray(zeros(NA+1), 0:NA)
global V_new = zeros(NA)


# variables to communicate with function
global ia_com

## Value function iteration and numerical minimization
function utility(x_in)
    # calculate consumption
    cons = max(a[ia_com] - x_in, 1e-10)

    # calculate future utility
    vplus = max(spline_eval(x_in), 1e-10)^egam/egam

    # get utility function

    return -(cons^egam/egam + beta*vplus)
end



# initialize a, c and value function
#global a = OffsetArray(a0.*Array(0:NA)./(NA+1), 0:NA)
#global c = OffsetArray(a/2.0, 0:NA)
#global a = a0.*Array(0:NA-1)./NA
global a = Array(LinRange(0.0, a0, NA))
global c = a/2.0

global spline_eval = CubicSpline(a,c)


# iterate until value function converges
for iter in 2:itermax
    # set a = 0 manually
    c[1] = 0.0
    V_new[1] = -1e10

    # calculate optimal decision for every gridpoint
    for ia in 1:NA

        # initialize starting value and communicate resources
        x_in = a[ia] - c[ia]
        ia_com = ia

        resultado = optimize(utility,0.0, a[ia], Brent())
        # get optimal consumption and value function
        c[ia] = a[ia] - resultado.minimizer
        V_new[ia] = -resultado.minimum

    end

    # Interpolate coefficients
    spline_eval = CubicSpline(a, (egam*V_new).^(1.0/egam)) 

    # get convergence level
    con_lev = maximum(@. abs(V_new - V )/max(abs(V), 1e-10))

    println(string(iter)*"   "*string(con_lev))

    #  check for convergence
    if con_lev < sig
        println("Convergio")
        break
    end

    V = copy(V_new)
end

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
plot!(1:TT,  c_t, label = "analytical")
