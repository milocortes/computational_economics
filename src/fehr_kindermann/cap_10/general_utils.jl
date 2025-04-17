#=#############################################################################
# SUBROUTINE discretize_AR
#
# Discretizes an AR(1) process of the form z_j = \rho*z_{j-1} + eps using
#     the Rouwenhorst method.
#
# REFERENCE: Kopecky, K.A., Suen, R.M.H., Finite state Markov-chain 
#            approximations to highly persistent processes, Review of Economic
#            Dynamics, Vol. 13, No. 3, 2010, 701-714.
=##############################################################################

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ)

    state_values, p = _rouwenhorst(p, p, m, ψ, N)
    #return p, state_values

    sigma_eta = σ /(1.0-ρ^2)
    psi_val = sqrt(N-1)*sqrt(sigma_eta)

    w_est = zeros(N)
    w_est = 1.0/float(N)
    for in in 1:10000
        w_est = *(transpose(p), w_est)
    end
 
    # [-psi_val + 2.0*psi_val*float(i-1)/float(n-1) for i in 1:N]
    return p , [-psi_val + (2.0*psi_val* ((i-1)/(N-1))) for i in 1:N], w_est
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)
        θN = p    *[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
             (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
             q    *[zeros(1, n); zeros(n-1, 1) θ_nm1] +
             (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

        θN[2:end-1, :] ./= 2

        return range(m-Δ, stop=m+Δ, length=n), θN
    end
end


###########################################################################
# FUNCTION get_tomorrow
#
# Calculates value of function that should be integrated for pis.
###########################################################################
function get_tomorrow(pi)

    ###### INPUT/OUTPUT VARIABLES #########################################

    # transition probabilities
    #real*8, intent(in) :: pi(:)

    # tomorrows shock
    #integer :: get_tomorrow


    ###### OTHER VARIABLES ################################################

    #real*8 :: rand
    #integer :: i1


    ###### ROUTINE CODE ###################################################

    # get random number
    #call random_number(rand)

    random_generated = rand()
    # get tomorrows value
    for i1 in 1:size(pi, 1)-1

        if(random_generated <= sum(pi[1:i1]))
            return i1
        end
    end

    # else choose last value
    return size(pi, 1)-1

end 


###############################################################################
# SUBROUTINE simulate_AR
#
# Simulates a discrete AR(1) process.
###############################################################################
function simulate_AR(pi, shocks)

        
    ###### INPUT/OUTPUT VARIABLES #############################################
    
    # transition matrix
    #real*8, intent(in) :: pi(:, :)
    
    # simulated schocks
    #integer, intent(out) :: shocks(:)
    
    # should the random seed be initialized at a fixed values
    #logical, optional :: fixed
    
    
    ###### OTHER VARIABLES ####################################################
    
    #integer :: T, n, j
    
    
    ###### ROUTINE CODE #######################################################
    
    # assert size equality and get number of simulated schocks
    #n = assert_eq(size(pi,1), size(pi,2), 'tauchen')
    n = size(pi,1)
    T = size(shocks,1)
    
    #= initialize the random seed
    if(tbox_seed)then
        if(present(fixed))then
            call init_random_seed(fixed)
        else
            call init_random_seed()
        endif
        tbox_seed = .false.
    endif
    =#

    # get first entry
    shocks[1] = round(Int, n ÷ (2+1))  
    
    # now calculate other shocks
    for j in 2:T
        shocks[j] = get_tomorrow( pi[round(Int,shocks[j-1]) , :] )
    end


    return shocks
end 

###############################################################################
# SUBROUTINE grid_Cons_Grow
#
# Constructs a growing grid on [left, right].
###############################################################################
function grid_Cons_Grow(a, n, left, right, growth)
    ccall((:grid_Cons_Grow, "./fortran_wrappers/mod_julfort.so"),
                Cvoid,
                (Ptr{Float64},Ref{Int64},Ref{Float64}, Ref{Float64}, Ref{Float64}),
                a, n,  left, right,growth)
end

###############################################################################
# subroutine linint_Grow
#
# Calculates linear interpolant on a growing grid.
###############################################################################
function linint_Grow(x, left, right, growth, n, il, ir, phi)

    ccall((:linint_Grow, "./fortran_wrappers/mod_julfort.so"),
            Cvoid,
            (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64} ),
            x, left, right, growth, n, il, ir, phi)

    return il[1], ir[1], phi[1]
end