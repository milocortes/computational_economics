###############################################################################
# PROGRAM SOLG_TR
#
# ## The stochastic OLG model with transitional dynamics
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################

include("utils_low_high_skills.jl")

using OffsetArrays
using Plots

# number of transition periods
global TT = 40

# Skills
global SS = 2

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
for param = [:r, :rn, :w, :wn, :p, :KK, :AA, :BB, :LL, :HH, :YY, :CC, :II, :GG, :INC]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

# government variables
for param = [:tauc, :tauw, :taur, :taup, :kappa, :tax, :PP]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

global gy
global by

global pen = OffsetArray(zeros(JJ, TT+1), 1:JJ,0:TT)

global taxrev = OffsetArray(zeros(4, TT+1), 1:4,0:TT) 

# LSRA variables
for param = [:BA, :SV]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

global lsra_comp
global lsra_all
global Lstar
global lsra_on

# cohort aggregate variables
for param = [:c_coh, :l_coh, :y_coh, :a_coh, :v_coh, :VV_coh]
    @eval global $param = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)
end

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 3

# demographic and other model parameters
for param = [:m, :pop]
    @eval global $param = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)
end

global eff = zeros(JJ, SS)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)

for param = [:aplus, :c, :l, :phi, :VV, :v, :phi]
    @eval global $param = OffsetArray(zeros(JJ, SS, NA+1, NP, NS, TT+1), 1:JJ, 1:SS, 0:NA, 1:NP, 1:NS, 0:TT)
end


global FLC = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)

# numerical variables
global RHS = OffsetArray(zeros(JJ, SS, NA+1, NP, NS, TT+1), 1:JJ, 1:SS, 0:NA, 1:NP, 1:NS, 0:TT)
global EV = OffsetArray(zeros(JJ, SS, NA+1, NP, NS, TT+1), 1:JJ, 1:SS, 0:NA, 1:NP, 1:NS, 0:TT)

global ij_com
global ia_com
global ip_com
global is_com
global it_com
global cons_com
global lab_com
global DIFF = OffsetArray(zeros(TT+1), 0:TT)


global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

# calculate initial equilibrium
get_SteadyState()

# set reform parameter (adjsust accordingly for Figure 11.7)
kappa[1:TT] .= 0.0

# calculate transition path without lsra
lsra_on = false
get_transition()