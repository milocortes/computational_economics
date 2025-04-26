###############################################################################
# PROGRAM SOLG_LR
#
# ## Long-run equilibria in the stochastic OLG model
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################
include("sol_prog_11_01_utils.jl")

using OffsetArrays
using Plots 

# number of years the household lives
global JJ = 12

# number of years the household retires
global JR = 10

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 5

# number of points on the asset grid
global NA = 100

# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0/gamma
global nu    = 0.335
global beta  = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps   = 0.05
global rho         = 0.98

# production parameters
global alpha = 0.36
global delta = 1.0-(1.0-0.0823)^5
global Omega = 1.60

# size of the asset grid
global a_l    = 0.0
global a_u    = 35.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+0.01)^5-1e0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 50


# macroeconomic variables
for param =[:r, :rn ,:w ,:wn ,:p ,:KK ,:AA ,:BB ,:LL ,:HH ,:YY ,:CC ,:II ,:GG ,:INC]
    @eval global $param
end

# government variables
for param =[:tauc, :tauw, :taur, :taup, :kappa, :gy, :by, :PP, :tax, :reform_on]
    @eval global $param
end

global pen = zeros(JJ)
global taxrev = zeros(4)

# cohort aggregate variables
for param =[:c_coh, :y_coh, :l_coh, :a_coh, :v_coh]
    @eval global $param = zeros(JJ)
end

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)
global is_initial = 3

# demographic and other model parameters
global m = zeros(JJ)
global eff = zeros(JJ)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)

for param = [:aplus, :c, :l, :phi, :V]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
end


# numerical variables
for param = [:RHS, :EV]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS)
end

for param = [:ij_com, :ia_com, :ip_com, :is_com, :it_com, :cons_com, :lab_com, :DIFF, :INC_init]
    @eval global $param
end

global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

initialize()

# calculate initial equilibrium
reform_on = false 

get_SteadyState()
