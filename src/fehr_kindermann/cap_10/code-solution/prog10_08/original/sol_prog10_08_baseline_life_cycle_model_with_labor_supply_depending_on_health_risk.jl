###############################################################################
# PROGRAM LaborHealth
#
# ## The baseline life cycle model with labor supply depending on health risk
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################
include("sol_prog10_08_utils.jl")

using OffsetArrays
using Plots 

# number of years the household lives
global JJ = 80

# year the household retires
global JR = 45

# number of points on the asset grid
global NA = 200

# number of health states
global NM = 1

# number of transitory shock process values
global NS = 7

# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0/gamma
global nu = 0.335
global beta = 0.98

# health cost parameters
global phi_l = 0.0
global chi = 0.80
global varrho = 0.20

# household risk process
global sigma_theta = 0.242
global sigma_eps   = 0.022
global rho         = 0.985

# size of the asset grid
global a_l    = 0.0
global a_u    = 200.0
global a_grow = 0.05

# probability of bad health
global dist_m = OffsetArray(zeros(NM+1), 0:NM)
global pi_m = OffsetArray(zeros(JJ, NM+1, NM+1), 1:JJ, 0:NM, 0:NM)

# the transitory shock process
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 4

# net prices
global r
global w

# transfer payments (old-age), survival probabilities, productivity
global pen = zeros(JJ)
global psi = OffsetArray(zeros(JJ+1, NM +1), 1:JJ+1, 0:NM)
global eff = zeros(JJ)
global varrho_m = OffsetArray(zeros(NM+1), 0:NM)

# individual variables

global a = OffsetArray(zeros(NA+1), 0:NA)

for param = [:aplus, :c, :l, :phi, :V]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NM+1, NS), 1:JJ, 0:NA, 0:NM, 1:NS)
end

# cohort aggregate variables

for param = [:c_coh, :l_coh, :y_coh, :a_coh, :v_coh, :cv_c, :cv_y, :frac_phi]
    @eval global $param = OffsetArray(zeros(JJ, NM+2),  1:JJ, 0:NM+1)
end


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NM+1, NS), 1:JJ, 0:NA, 0:NM, 1:NS)
global EV = OffsetArray(zeros(JJ, NA+1, NM+1, NS), 1:JJ, 0:NA, 0:NM, 1:NS)

global ij_com
global ia_com
global im_com
global is_com
global cons_com
global lab_com

global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

# initialize remaining variables
@elapsed initialize()

# solve the household problem
@elapsed solve_household()

# calculate the distribution of households over state space
@elapsed get_distribution()

# aggregate individual decisions
@elapsed aggregation()

ages = 20 .+ (1:JJ)

plot(ages, c_coh[:, 0], ylabel = "Consumption", xlabel = "Age j", label="Good Health")
plot!(ages, c_coh[:, 1], label="Bad Health")
plot!(ages, c_coh[:, 2], label="Average")


plot(ages, l_coh[:, 0], ylabel = "Working time", xlabel = "Age j", label="Good Health")
plot!(ages, l_coh[:, 1], label="Bad Health")
plot!(ages, l_coh[:, 2], label="Average")

plot(ages, y_coh[:, 0], ylabel = "Earnings", xlabel = "Age j", label="Good Health")
plot!(ages, y_coh[:, 1], label="Bad Health")
plot!(ages, y_coh[:, 2], label="Average")
