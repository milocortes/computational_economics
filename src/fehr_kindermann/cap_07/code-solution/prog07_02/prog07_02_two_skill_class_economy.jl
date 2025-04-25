###############################################################################
# PROGRAM SKILL_OLG
#
# ##  Two skill class economy
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################
include("prog07_02_utils.jl")

using OffsetArrays
using Plots 

# model parameters
global TT = 24    
global JJ = 3     
global JR = 3     
global SS = 2
global  gamma = 0.5
global  egam = 1.0 - 1.0/gamma
global  beta = 0.9
global  nu = 1.5
global  rho = 0.6
global  erho = 1.0 - 1.0/rho
global  alpha = 0.3
global  delta = 0.0
global  tol = 1e-6
global  damp = 0.2
global itermax = 200


# model variables
## Macro variables
for param = [:w , :r, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :gy, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, :Tpen, :TXR, :tax, :eps]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
end

## Cohort variables
for param = [:mu, :wn, :m, :a, :c, :util, :l, :pen]
    @eval global $param = OffsetArray(zeros(JJ, SS, TT+1), 1:JJ, 1:SS, 0:TT)
end

global h = zeros(JJ, SS)
global v = OffsetArray(zeros(SS, TT+JJ-2+1), 1:SS,-JJ+2:TT)

global lsra_on
global smopec

