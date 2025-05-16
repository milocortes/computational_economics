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

include("utils_sol_prog_11_02_SAM.jl")

using OffsetArrays
using Plots
using DataFrames
using CSV 

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
global sig     = 1e-5
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
    @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
end

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi = zeros(NS, NS)
global eta = zeros(NS)

global is_initial = 3

# demographic and other model parameters
for param = [:m, :pop]
    @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT)
end

global eff = zeros(JJ)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)

for param = [:aplus, :c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT)
end


global FLC = OffsetArray(zeros(JJ, TT+1), 1:JJ,0:TT)

# numerical variables
global RHS = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT)
global EV = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT)

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


# Compute social accounts
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
    PEN = [pen[JR, 0], kappa[0]],
    PP = [PP[0], PP[0]/YY[0]*100]
);

capital_market
labour_market
good_market
gov_accounts
pension_system

data = ["Producción" NaN NaN NaN NaN CC[0] NaN GG[0] II[0]  NaN;
"Salarios" w[0]*LL[0]  NaN NaN NaN NaN NaN NaN NaN NaN;
"Ganancias" (r[0]+delta)KK[0] NaN NaN NaN NaN NaN NaN NaN NaN;
"Pensiones" NaN NaN NaN NaN w[0]*LL[0]*taup[0]  NaN NaN NaN NaN;
"Hogares" NaN w[0]*LL[0]*(1-tauw[0]) (r[0]+delta)*KK[0]*(1-taur[0]) sum(pen[:,0].*m[:,0]) NaN NaN BB[0]*r[0] NaN NaN;
"Gobierno (Ingresos)" NaN w[0]*LL[0]*taur[0] (r[0]+delta)*KK[0]*taur[0] NaN tauc[0]*CC[0] NaN NaN NaN NaN;
"Gobierno (Gasto)" NaN NaN NaN NaN NaN GG[0]+(1+r[0])*BB[0] - (1+n_p)*BB[0]+BB[0]*r[0]  NaN  NaN  NaN;
"Ahorro-Inversión" NaN NaN NaN NaN r[0]*AA[0]+BB[0]*r[0] ( (1+r[0])*BB[0]*n_p) (r[0] - n_p)*BB[0] NaN NaN;
"Total" NaN NaN NaN NaN NaN NaN NaN NaN NaN]


SAM = DataFrame(
     data,
    ["", "Producción", "Salarios", "Ganancias", "Pensiones", "Hogares", "Gobierno (Ingresos)", "Gobierno (Gasto)", "Ahorro-Inversión", "Total"]

)

CSV.write("outputs/sam_olg.csv", SAM)

# set reform parameter (adjsust accordingly for Figure 11.7)
kappa[1:TT] .= 0.0

# calculate transition path without lsra
lsra_on = false
get_transition()