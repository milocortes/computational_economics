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
include("utils_persistence_shocks.jl")

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

# number of persistent shock process values
global NM = 1

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

# health cost parameters
global phi_l = 0.0
global chi = 0.80
global varrho = 0.20

# production parameters
global alpha = 0.36
global delta = 1.0-(1.0-0.0823)^5
global Omega2 = 1.10

# size of the asset grid
global a_l    = 0.0
global a_u    = 50.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+0.02)^5-1.0

# simulation parameters
global damp    = 0.20
global sig     = 1e-4
global itermax = 70

# counter variables
#global iter

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

# LSRA variables
global BA = OffsetArray(zeros(TT+1), 0:TT)
global SV = OffsetArray(zeros(TT+1), 0:TT)
global lsra_comp
global lsra_all
global Vstar
global lsra_on

# cohort aggregate variables
for param = [:c_coh, :l_coh, :y_coh, :a_coh, :v_coh, :VV_coh, :frac_phi]
    @eval global $param = OffsetArray(zeros(JJ, NM +1, TT+1), 1:JJ, 0:NM, 0:TT )
end

for param = [:GAM, :beq]
    @eval global $param = OffsetArray(zeros(JJ, NM+1, TT+1), 1:JJ, 0:NM, 0:TT )
end

global omega = zeros(JJ)
global beq_coh = OffsetArray(zeros(JJ, NM+1, TT+1) , 1:JJ, 0:NM, 0:TT)
global psi = OffsetArray(zeros(JJ+1, NM+1, TT+1), 1:JJ+1, 0:NM, 0:TT)
global varrho_m = OffsetArray(zeros(NM+1), 0:NM)

# the shock process
global dist_theta = OffsetArray(zeros(NM+1), 0:NM)
global theta = OffsetArray(zeros(NM+1), 0:NM)
global pi = zeros(NS, NS)
global eta = zeros(NS)
global is_initial = 4

# demographic and other model parameters
global m = OffsetArray(zeros(JJ, NM +1, TT+1), 1:JJ, 0:NM, 0:TT)
#global eff = OffsetArray(zeros(JJ, NM+1), 1:JJ, 0:NM)
global eff = zeros(JJ)

# probability of bad health
global dist_m = OffsetArray(zeros(NM+1), 0:NM)
global pi_m = OffsetArray(zeros(JJ, NM+1, NM+1), 1:JJ, 0:NM, 0:NM)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)
global FLC = OffsetArray(zeros(JJ, NM +1, TT+1), 1:JJ, 0:NM, 0:TT)

for param = [:aplus, :c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NM +1, NS, TT+1), 1:JJ, 0:NA, 0:NM, 1:NS, 0:TT)
end


# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NM +1, NS, TT+1), 1:JJ, 0:NA, 0:NM, 1:NS, 0:TT)
global EV = OffsetArray(zeros(JJ, NA+1, NM +1, NS, TT+1), 1:JJ, 0:NA, 0:NM, 1:NS, 0:TT)
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


capital_market = DataFrame( 
          etiqueta = ["valor", "(in %)"],
          K = [KK[0], KK[0]/YY[0]*500],
          A = [AA[0], AA[0]/YY[0]*500],  
          B = [BB[0], BB[0]/YY[0]*500], 
          BA = [BA[0], BA[0]/YY[0]*500],   
          r = [r[0], ""],
          pa = [((1.0+r[0])^(1.0/5.0)-1.0)*100.0,  ""]
);

labour_market = DataFrame( 
    etiqueta = ["valor"],
    L = [LL[0]],
    HH = [HH[0]*100],
    INC = [INC[0]],
    w = [w[0]]
);

good_market = DataFrame( 
          etiqueta = ["valor", "(in %)"],
          Y = [YY[0], YY[0]/YY[0]*100],
          C = [CC[0], CC[0]/YY[0]*100],  
          I = [II[0], II[0]/YY[0]*100], 
          G = [GG[0], GG[0]/YY[0]*100]
); 

gov_accounts = DataFrame(
    etiqueta = ["valor", "(in %)", "(rate)"],
    TAUC = [taxrev[1,0], taxrev[1,0]/YY[0]*100, tauc[0]*100],   
    TAUW = [taxrev[2,0], taxrev[2,0]/YY[0]*100, tauw[0]*100],
    TAUR = [taxrev[3,0], taxrev[3,0]/YY[0]*100, taur[0]*100],
    TOTAL = [taxrev[4,0], taxrev[4,0]/YY[0]*100, ""], 
    G = [GG[0], GG[0]/YY[0]*100, ""],
    B = [BB[0], (BB[0]*5.0)/YY[0]*100, ""]
);

pension_system = DataFrame(
    etiqueta = ["valor", "(in %)"],
    TAUP = [taup[0]*w[0]*LL[0], taup[0]*100], 
    PEN = [pen[JR, 0] , kappa[0]],
    PP = [PP[0], PP[0]/YY[0]*100]
);

capital_market
labour_market
good_market
gov_accounts
pension_system


# set reform parameters
kappa[1:TT] .= 0.7

# calculate transition path without lsra
lsra_on = false
get_transition()



# The Long-Run Eﬀect
## Long-run effects of the consumption tax reform over the life cycle of the households.
### Private Consumption
plot([i for i in 20:5:75],  parent(mean(c_coh[1:12,:, 0], dims=2))  , title = "Long-run Effects on Consumption", xlabel = "Year", label = "Consumption - Pre-Reforma")
plot!([i for i in 20:5:75],   parent(mean(c_coh[1:12,:, 140], dims=2))  , label = "Consumption- Post-Reforma")

### Working Hours
plot([i for i in 20:5:75],  parent(mean(l_coh[1:12,:, 0], dims=2)), title = "Long-run Effects on Hours Worked", xlabel = "Year", label = "Pre-Reforma")
plot!([i for i in 20:5:75],  mean(l_coh[1:12,:, 140], dims=2), label = "Post-Reforma")

### Earnings
plot([i for i in 20:5:75],  parent(mean(w[0].*y_coh[1:12,:, 0], dims=2)), title = "Average life-cycle", label = "Earnings - Pre-Reforma")
plot!([i for i in 20:5:75],  mean(w[0].*y_coh[1:12,:, TT], dims=2), label = "Earnings - Post-Reforma")

### Private Wealth
plot([i for i in 20:5:75],  parent(mean(a_coh[1:12,:, 0], dims=2)), title = "Average life-cycle", label = "Wealth - Pre-Reforma")
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


# calculate transition path with lsra
lsra_on = true
get_transition()
