###############################################################################
# PROGRAM SOLG_TR_PEN
#
# ## Implicit tax rates in the pension system
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################

using OffsetArrays

# number of transition periods
global TT = 40

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

# number of points on the earnings points grid for retirement
global NR = 10

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

# size of the earnings point grid
global ep_l    = 0.0
global ep_u    = 7.0
global ep_grow = 0.02

# demographic parameters
global n_p   = (1.0+0.01)^5-1.0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 50

# counter variables
#integer :: iter

# macroeconomic variables
for param = [:r, :rn, :w, :wn, :p, :KK, :AA, :BB, :LL, :HH, :YY, :CC, :II, :GG]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

# government variables
global gy
global by

for param = [:tauc , :tauw, :taur, :taup, :kappa, :lambda, :PP, :tax]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

global tau_impl = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
global taxrev = OffsetArray(zeros(4,TT+1), 1:4,0:TT)


# LSRA variables
for param = [:BA, :SV]
    @eval global $param = OffsetArray(zeros(TT+1),0:TT)    
end

global lsra_comp
global lsra_all
global Lstar
global lsra_on

# cohort aggregate variables
for param = [:c_coh, :l_coh, :y_coh, :pen, :v_coh, :VV_coh]
    @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
end

global a_coh = OffsetArray(zeros(JJ+1, TT+1), 1:JJ+1, 0:TT)



# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 3

# demographic and other model parameters
global m = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
global pop = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
global workpop = OffsetArray(zeros(TT+1), 0:TT)
global INC = OffsetArray(zeros(TT+1), 0:TT)

global eff = zeros(JJ)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)
global ep = OffsetArray(zeros(NR+1), 0:NR)

for param = [:aplus, :epplus, :c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NR+1, NP, NS, TT+1), 1:JJ, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)
end

global penp = OffsetArray(zeros(JJ, TT+1, NR+1), 1:JJ, 0:TT, 0:NR)
global FLC = OffsetArray(zeros(JJ, TT+1), 1:JJ,0:TT)


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NR+1, NP, NS, TT+1), 1:JJ, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)
global EV = OffsetArray(zeros(JJ, NA+1, NR+1, NP, NS, TT+1), 1:JJ, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)

for param = [:ij_com, :ia_com, :ir_com, :ip_com, :is_com, :it_com, :cons_com, :lab_com, :epplus_com]
    @eval global $param
end

global DIFF = OffsetArray(zeros(TT+1), 0:TT)


## Auxiliares para grid de ahorro
global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

## Auxiliares para grid de earning points
global ial_ep = Array{Int64}(undef, 1)
global iar_ep = Array{Int64}(undef, 1)
global varphi_ep = zeros(1)
