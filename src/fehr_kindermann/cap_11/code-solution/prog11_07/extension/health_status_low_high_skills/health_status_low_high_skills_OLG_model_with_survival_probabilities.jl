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
include("health_status_utils.jl")

using OffsetArrays
using Plots 

# number of transition periods
global TT = 40

# number of years the household lives-check
global JJ = 16

# year the household retires-check
global JR = 10

# number of points on the asset grid-check
global NA = 200

# number of health states-check
global NM = 1

# number of transitory shock process values-check
global NS = 7

# number of persistent shock process values
global NP = 2

# number of health shock process values
global NH = 5


# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0/gamma
global beta = 0.98

# household risk process
global sigma_theta = 0.242
global sigma_eps   = 0.022
global rho         = 0.985

# the health process
global delta      = 0.1
global chi        = 1.0
global sigma_zeta = 0.2
global c_floor    = 0.15

# production parameters
global alpha = 0.36
global delta = 1.0-(1.0-0.0823)^5
global Omega2 = 1.60

# size of the asset grid
global a_l    = 0.0
global a_u    = 200.0
global a_grow = 0.05

# the permanent shock process
global dist_theta = zeros(NP);
global theta = zeros(NP);

# the health cost shock
global  dist_zeta = zeros(NH);
global zeta = zeros(NH);
global  hc = zeros(JJ, NH);
global varrho = zeros(NP);


global varrho_m = OffsetArray(zeros(NP, NM+1), 1:NP, 0:NM);


# demographic parameters
global n_p   = (1.0+0.02)^5-1.0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 70

# probability of bad health
global dist_m = OffsetArray(zeros(NM+1), 0:NM)
global pi_m = OffsetArray(zeros(JJ, NP, NM+1, NM+1), 1:JJ, 1:NP, 0:NM, 0:NM);

# the transitory shock process
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 4

# macroeconomic variables
for param = [:r, :rn, :w, :wn, :p, :KK, :AA, :BB, :LL, :HH, :YY, :CC, :II, :GG, :INC, :BQ]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end


# government variables
for param = [:tauc, :tauw, :taur, :taup, :kappa, :PP, :tax]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

global taxrev = OffsetArray(zeros(4, TT+1), 1:4,0:TT)
global gy
global by

global pen = OffsetArray(zeros(JJ, TT+1), 1:JJ,0:TT)


# transfer payments (old-age), survival probabilities, productivity
#global pen = zeros(JJ)
global psi = OffsetArray(zeros(JJ+1, NM +1, TT +1), 1:JJ+1, 0:NM, 0:TT)
global eff = zeros(JJ)
global a_bor = zeros(JJ, NP);
global varrho_m = OffsetArray(zeros(NM+1), 0:NM)

# demographic and other model parameters
global m = OffsetArray(zeros(JJ, NM+1, TT+1), 1:JJ, 0:NM, 0:TT)

# individual variables

global a = OffsetArray(zeros(NA+1), 0:NA)

for param = [:aplus, :c, :l, :phi, :V, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS, TT+1), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS, 0:TT)
end

# cohort aggregate variables

for param = [:c_coh, :l_coh, :y_coh, :a_coh, :v_coh, :cv_c, :cv_y, :frac_phi, :GAM, :beq]
    @eval global $param = OffsetArray(zeros(JJ, NP+1, NM+2, TT+1),  1:JJ, 1:NP+1, 0:NM+1, 0:TT)
end

global omega = zeros(JJ)
global beq_coh = OffsetArray(zeros(JJ, NM+1, TT+1) , 1:JJ, 0:NM, 0:TT)
global beq_generation = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)

# OUT OF THE POKET
global k = zeros(JJ);


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS);
global EV = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS);
global ij_com
global ia_com
global ip_com
global im_com
global is_com
global ih_com
global cons_com

global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

# calculate initial equilibrium
get_SteadyState()

# set reform parameters
kappa[1:TT] .= 0.5

# calculate transition path without lsra
lsra_on = false
get_transition()
