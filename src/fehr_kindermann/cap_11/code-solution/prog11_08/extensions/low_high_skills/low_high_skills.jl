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
include("utils_low_high_skills.jl")

using OffsetArrays
using DataFrames
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

# number of points on the earnings points grid for retirement
global NR = 10

# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0 / gamma
global nu = 0.335
global beta = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps = 0.05
global rho = 0.98

# production parameters
global alpha = 0.36
global delta = 1.0 - (1.0 - 0.0823)^5
global Omega = 1.60

# size of the asset grid
global a_l = 0.0
global a_u = 55.0
global a_grow = 0.05

# size of the earnings point grid
global ep_l = 0.0
global ep_u = 10.0
global ep_grow = 0.02

# demographic parameters
global n_p = (1.0 + 0.01)^5 - 1.0

# simulation parameters
global damp = 0.30
global sig = 1e-4
global itermax = 50

# counter variables
#integer :: iter

# macroeconomic variables
for param = [:r, :rn, :w, :wn, :p, :KK, :AA, :BB, :LL, :HH, :YY, :CC, :II, :GG]
    @eval global $param = OffsetArray(zeros(TT + 1), 0:TT)
end

# government variables
global gy
global by

for param = [:tauc, :tauw, :taur, :taup, :kappa, :lambda, :PP, :tax]
    @eval global $param = OffsetArray(zeros(TT + 1), 0:TT)
end

global tau_impl = OffsetArray(zeros(JJ, TT + 1), 1:JJ, 0:TT)
global taxrev = OffsetArray(zeros(4, TT + 1), 1:4, 0:TT)


# LSRA variables
for param = [:BA, :SV]
    @eval global $param = OffsetArray(zeros(TT + 1), 0:TT)
end

global lsra_comp
global lsra_all
global Lstar
global lsra_on

# cohort aggregate variables
for param = [:c_coh, :l_coh, :y_coh, :a_coh, :pen, :v_coh, :VV_coh]
    @eval global $param = OffsetArray(zeros(JJ, SS, TT + 1), 1:JJ, 1:SS, 0:TT)
end


# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 3

# demographic and other model parameters
global m = OffsetArray(zeros(JJ, SS, TT + 1), 1:JJ, 1:SS, 0:TT)
global pop = OffsetArray(zeros(JJ, SS, TT + 1), 1:JJ, 1:SS, 0:TT)
global workpop = OffsetArray(zeros(TT + 1), 0:TT)
global INC = OffsetArray(zeros(TT + 1), 0:TT)

global eff = zeros(JJ, SS)

# individual variables
global a = OffsetArray(zeros(NA + 1), 0:NA)
global ep = OffsetArray(zeros(NR + 1), 0:NR)

for param = [:aplus, :epplus, :c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, SS, NA + 1, NR + 1, NP, NS, TT + 1), 1:JJ, 1:SS, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)
end

global penp = OffsetArray(zeros(JJ, SS, TT + 1, NR + 1), 1:JJ, 1:SS, 0:TT, 0:NR)
global FLC = OffsetArray(zeros(JJ, SS, TT + 1), 1:JJ, 1:SS, 0:TT)


# numerical variables
global RHS = OffsetArray(zeros(JJ, SS, NA + 1, NR + 1, NP, NS, TT + 1), 1:JJ, 1:SS, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)
global EV = OffsetArray(zeros(JJ, SS, NA + 1, NR + 1, NP, NS, TT + 1), 1:JJ, 1:SS, 0:NA, 0:NR, 1:NP, 1:NS, 0:TT)

for param = [:ij_com, :ia_com, :ir_com, :ip_com, :is_com, :it_com, :cons_com, :lab_com, :epplus_com]
    @eval global $param
end

global DIFF = OffsetArray(zeros(TT + 1), 0:TT)

## Auxiliares para grid de ahorro
global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

## Auxiliares para grid de earning points
global ial_ep = Array{Int64}(undef, 1)
global iar_ep = Array{Int64}(undef, 1)
global varphi_ep = zeros(1)

## Files
global file_output
global file_summary

# calculate initial equilibrium
get_SteadyState()

# Compute social accounts
capital_market = DataFrame(
    etiqueta=["valor", "(in %)"],
    K=[KK[0], KK[0] / YY[0] * 500],
    A=[AA[0], AA[0] / YY[0] * 500],
    B=[BB[0], BB[0] / YY[0] * 500],
    BA=[BA[0], BA[0] / YY[0] * 500],
    r=[r[0], ""],
    pa=[((1.0 + r[0])^(1.0 / 5.0) - 1.0) * 100.0, ""]
);

labour_market = DataFrame(
    etiqueta=["valor"],
    L=[LL[0]],
    HH=[HH[0] * 100],
    INC=[INC[0]],
    w=[w[0]]
);
good_market = DataFrame(
    etiqueta=["valor", "(in %)"],
    Y=[YY[0], YY[0] / YY[0] * 100],
    C=[CC[0], CC[0] / YY[0] * 100],
    I=[II[0], II[0] / YY[0] * 100],
    G=[GG[0], GG[0] / YY[0] * 100]
);

gov_accounts = DataFrame(
    etiqueta=["valor", "(in %)", "(rate)"],
    TAUC=[taxrev[1, 0], taxrev[1, 0] / YY[0] * 100, tauc[0] * 100],
    TAUW=[taxrev[2, 0], taxrev[2, 0] / YY[0] * 100, tauw[0] * 100],
    TAUR=[taxrev[3, 0], taxrev[3, 0] / YY[0] * 100, taur[0] * 100],
    TOTAL=[taxrev[4, 0], taxrev[4, 0] / YY[0] * 100, ""],
    G=[GG[0], GG[0] / YY[0] * 100, ""],
    B=[BB[0], (BB[0] * 5.0) / YY[0] * 100, ""]
);

pension_system = DataFrame(
    etiqueta=["valor", "(in %)"],
    TAUP=[taup[0] * w[0] * LL[0], taup[0] * 100],
    #PEN=[pen[JR, 0], kappa[0]],
    PP=[PP[0], PP[0] / YY[0] * 100]
);

# set reform parameters (adjust accordingly for Figure 11.8)
#lambda(0:TT) = 0.99d0
#kappa[1:TT] .= 0.19
kappa[1:TT] .= 0.7

# calculate transition path without lsra
lsra_on = false
get_transition()

# The Long-Run Eﬀect
## Long-run effects of the consumption tax reform over the life cycle of the households.
### Private Consumption
plot([i for i in 20:5:75],  mean(c_coh[1:12,:, 0], dims=2)  , title = "Long-run Effects on Consumption", xlabel = "Year", label = "Consumption - Pre-Reforma")
plot!([i for i in 20:5:75],   mean(c_coh[1:12,:, TT], dims=2)  , label = "Consumption- Post-Reforma")

### Working Hours
plot([i for i in 20:5:75],  mean(l_coh[1:12,:, 0], dims=2), title = "Long-run Effects on Hours Worked", xlabel = "Year", label = "Pre-Reforma")
plot!([i for i in 20:5:75],  mean(l_coh[1:12,:, TT], dims=2), label = "Post-Reforma")

### Earnings
plot([i for i in 20:5:75],  mean(w[0].*y_coh[1:12,:, 0], dims=2), title = "Average life-cycle", label = "Earnings - Pre-Reforma")
plot!([i for i in 20:5:75],  mean(w[0].*y_coh[1:12,:, TT], dims=2), label = "Earnings - Post-Reforma")

### Private Wealth
plot([i for i in 20:5:75],  mean(a_coh[1:12,:, 0], dims=2), title = "Average life-cycle", label = "Wealth - Pre-Reforma")
plot!([i for i in 20:5:75],  mean(a_coh[1:12,:, TT], dims=2), label = "Wealth Worked - Post-Reforma")



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

df_transition_path = first(df_transition_path, TT+1)

image_transition_path = plot( plot(df_transition_path.labour_supply, title = "Labour supply"), 
plot(df_transition_path.gdp, title = "Gross domestic product"),  
plot(df_transition_path.consumption, title = "Private consumption"),
plot(df_transition_path.interest_rate, title = "Interest Rate"),
plot(df_transition_path.wage, title = "Wage"),
plot(df_transition_path.pensiones, title = "Pensions (%GDP)"),

layout = (3, 2), ylabel = "%Ch from baseline", yguidefontsize = 7);

image_transition_path

### Calculate transition path with lsra
lsra_on = true
get_transition()

# close files
close(file_output)
close(file_summary)

