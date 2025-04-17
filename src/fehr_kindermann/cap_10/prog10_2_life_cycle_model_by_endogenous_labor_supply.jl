include("prog10_2_utils.jl")

using OffsetArrays

# number of years the household retires
global JR = 45

# number of years the household lives
global JJ = 80

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 7

# number of points on the asset grid
global NA = 200

# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0/gamma
global nu    = 0.335
global beta  = 0.98

# household risk process
global sigma_theta = 0.242
global sigma_eps   = 0.022
global rho         = 0.985

# size of the asset grid
global a_l    = 0.0
global a_u    = 200.0
global a_grow = 0.05


# net prices
global r
global w

# transfer payments (old-age) and survival probabilities
global pen = zeros(JJ)
global psi = zeros(JJ+1)

# cohort aggregate variables
for param = [:c_coh, :y_coh, :l_coh, :h_coh, :a_coh, :v_coh, :cv_c, :cv_y, :cv_l, :cv_h, :corr_hl]
    @eval global $param = zeros(JJ)
end

# the permanent shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)

# the transitory shock process
global pi = zeros((NS, NS))
global eta = zeros(NS)

global is_initial = 4

# demographic and other model parameters
global eff = zeros(JJ)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)
global aplus = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
global c = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
global l = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
global phi = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)

global V = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS) 


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
global EV = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)

global ij_com
global ia_com
global ip_com
global is_com
global cons_com
global lab_com

global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)


# initialize remaining variables
initialize()

# solve the household problem
solve_household()

# calculate the distribution of households over state space
get_distribution()

# aggregate individual decisions
aggregation()

ages = 20 .+ (1:JJ)

plot(ages, c_coh, label="Consumption", xlabel="Age j", ylabel="Mean", title = "Extension of the life cycle model by\n endogenous labor supply (Julia!!)")
plot!(ages, y_coh+pen, label="Earnings")
plot!(ages, l_coh, label="Hours")

