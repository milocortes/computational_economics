###############################################################################
# PROGRAM SOLG_TR_SRV
#
# ## OLG model with survival probabilities
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################
include("utils_mex_skills.jl")

using OffsetArrays
using Plots 
using Statistics
using DataFrames

# number of transition periods
global TT = 140

# number of years the household lives
global JJ = 16

# number of years the household retires
global JR = 10

# skills classes
global SS = 2

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 5

# number of points on the asset grid
global NA = 100

# household preference parameters
global gamma = 0.18
global egam = 1.0 - 1.0/gamma
global nu    = 0.3890909090909091
global beta  = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps   = 0.05
global rho         = 0.98

# production parameters
global alpha = 0.622000128030777
global delta = 0.17767937625290053
global Omega2 = 1.45

# size of the asset grid
global a_l    = 0.0
global a_u    = 50.0
global a_grow = 0.05

# demographic parameters
global n_p   = 0.11167865656497411

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 70

# counter variables
global iter

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

global pen = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)

# LSRA variables
global BA = OffsetArray(zeros(TT+1), 0:TT)
global SV = OffsetArray(zeros(TT+1), 0:TT)
global lsra_comp
global lsra_all
global Vstar
global lsra_on

# cohort aggregate variables
for param = [:c_coh, :l_coh, :y_coh, :a_coh, :v_coh, :VV_coh, :GAM, :beq]
    @eval global $param = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT )
end


global omega = zeros(JJ)
global beq_coh = OffsetArray(zeros(JJ, SS, TT+1) , 1:JJ, 1:SS, 0:TT)
global psi = OffsetArray(zeros(JJ+1, TT+1), 1:JJ+1, 0:TT)

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)
global is_initial = 3

# demographic and other model parameters
global m = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)
global eff = zeros(JJ, SS)

# individual variables
global a = OffsetArray(zeros(NA+1, SS), 0:NA, 1:SS)
global FLC = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)

for param = [:aplus, :c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS, SS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 1:SS, 0:TT)
end


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NP, NS, SS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 1:SS, 0:TT)
global EV = OffsetArray(zeros(JJ, NA+1, NP, NS, SS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 1:SS, 0:TT)
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

# set reform parameters
#kappa[1:TT] .= 0.5
psi = min.(psi*1.05,1.0)
global sig     = 1e-4
global itermax = 10
# calculate transition path without lsra
lsra_on = false
get_transition()


# The Long-Run Eﬀect
## Long-run effects of the consumption tax reform over the life cycle of the households.
### Private Consumption
plot([i for i in 20:5:75],  mean(c_coh[1:12,:, 0], dims=2)  , title = "Long-run Effects on Consumption", xlabel = "Year", label = "Consumption - Pre-Reforma")
plot!([i for i in 20:5:75],   mean(c_coh[1:12,:, 140], dims=2)  , label = "Consumption- Post-Reforma")

### Working Hours
plot([i for i in 20:5:75],  mean(l_coh[1:12,:, 0], dims=2), title = "Long-run Effects on Hours Worked", xlabel = "Year", label = "Pre-Reforma")
plot!([i for i in 20:5:75],  mean(l_coh[1:12,:, 140], dims=2), label = "Post-Reforma")

### Earnings
plot([i for i in 20:5:75],  mean(w[0].*y_coh[1:12,:, 0], dims=2), title = "Average life-cycle", label = "Earnings - Pre-Reforma")
plot!([i for i in 20:5:75],  mean(w[0].*y_coh[1:12,:, TT], dims=2), label = "Earnings - Post-Reforma")

### Private Wealth
plot([i for i in 20:5:75],  mean(a_coh[1:12,:, 0], dims=2), title = "Average life-cycle", label = "Wealth - Pre-Reforma")
plot!([i for i in 20:5:75],  mean(a_coh[1:12,:, 140], dims=2), label = "Wealth Worked - Post-Reforma")


### The Transition Eﬀect
## Transition effect of the consumption tax reform on macroeconomic variables.
# Labor supply (in efficiency units)
(LL[0:TT]./LL[0].-1 )*100
# Gross domestic product
(YY[0:TT]./YY[0] .-1) .*100 
# Private consumption
(CC[0:TT]./CC[0] .-1) .*100 
# Interest rate
(r[0:TT]./r[0] .-1) .*100 
# Payroll taxes
((taup[0:TT].*w[0:TT].*LL[0:TT]./YY[0:TT])/(taup[0].*w[0].*LL[0]./YY[0]) .-1)*100
(((pen[JR, 1, 0:TT] + pen[JR, 2, 0:TT])./YY[0:TT])/((pen[JR, 1, 0] + pen[JR, 2, 0])./YY[0]).-1)*100
# Wage
(w[0:TT]./w[0] .-1) .*100 

df_transition_path = DataFrame(
    labour_supply = (LL[0:TT]./LL[0].-1 )*100,
    gdp = (YY[0:TT]./YY[0] .-1) .*100 ,
    consumption = (CC[0:TT]./CC[0] .-1) .*100 ,
    interest_rate = (r[0:TT]./r[0] .-1) .*100 ,
    wage = (w[0:TT]./w[0] .-1) .*100, 
    pensiones = (((pen[JR, 1, 0:TT] + pen[JR, 2, 0:TT])./YY[0:TT])/((pen[JR, 1, 0] + pen[JR, 2, 0])./YY[0]).-1)*100
)

df_transition_path = first(df_transition_path, 130)

image_transition_path = plot( plot(df_transition_path.labour_supply, title = "Labour supply"), 
plot(df_transition_path.gdp, title = "Gross domestic product"),  
plot(df_transition_path.consumption, title = "Private consumption"),
plot(df_transition_path.interest_rate, title = "Interest Rate"),
plot(df_transition_path.wage, title = "Wage"),
plot(df_transition_path.pensiones, title = "Pensions (%GDP)"),

layout = (3, 2), ylabel = "%Ch from baseline", yguidefontsize = 7);