###############################################################################
# PROGRAM SOLG_LR
#
# ## Long-run equilibria in the baseline life cycle model with health risk
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################


using OffsetArrays
using Plots 

# number of years the household lives
global JJ = 80

# year the household retires
global JR = 45

# number of points on the asset grid
global NA = 200

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 7

# number of health shock process values
global NH = 5

# maximum number of health states (plus 0)
global NM = 1

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

# size of the asset grid
global a_l    = 0.0
global a_u    = 600.0
global a_grow = 0.05

# the permanent shock process
global dist_theta = zeros(NP);
global theta = zeros(NP);

# the transitory shock process
global pi = zeros(NS, NS);
global eta = zeros(NS);
global is_initial = 4

# the health shock
global  dist_m = OffsetArray(zeros(NM+1), 0:NM);
global pi_m = OffsetArray(zeros(JJ, NP, NM+1, NM+1), 1:JJ, 1:NP, 0:NM, 0:NM);

# the health cost shock
global  dist_zeta = zeros(NH);
global zeta = zeros(NH);
global  hc = zeros(JJ, NH);
global varrho = zeros(NP);


global varrho_m = OffsetArray(zeros(NP, NM+1), 1:NP, 0:NM);

# net prices
global r
global w

# transfer payments (old-age), survival probabilities, productivity
global pen = zeros(JJ);
global psi = OffsetArray(zeros(JJ+1, NM+1), 1:JJ+1, 0:NM);
global eff = zeros(JJ);
global a_bor = zeros(JJ, NP);

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)
global aplus = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS, NH), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS, 1:NH);
global c = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS, NH), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS, 1:NH);
global phi = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS, NH), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS, 1:NH);
global V = OffsetArray(zeros(JJ, NA+1, NP, NM+1, NS, NH), 1:JJ, 0:NA, 1:NP, 0:NM, 1:NS, 1:NH);

# cohort aggregate variables
global c_coh = OffsetArray(zeros(JJ, NP+1, NM+2), 1:JJ, 1:NP+1, 0:NM+1);
global y_coh = OffsetArray(zeros(JJ, NP+1, NM+2), 1:JJ, 1:NP+1, 0:NM+1);
global a_coh = OffsetArray(zeros(JJ, NP+1, NM+2), 1:JJ, 1:NP+1, 0:NM+1);
global v_coh = OffsetArray(zeros(JJ, NP+1, NM+2), 1:JJ, 1:NP+1, 0:NM+1);
global cv_c = zeros(JJ);
global cv_y = zeros(JJ);
global frac_phi = OffsetArray(zeros(JJ, NP, NM+2), 1:JJ, 1:NP, 0:NM+1);

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
