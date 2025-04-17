include("general_utils.jl")
using Roots

# the first order condition
function foc(x_in)
    global ij_com 
    global ia_com 
    global is_com 
    global lab_com 
    global cons_com 

    # calculate tomorrows assets
    a_plus  = x_in

    # calculate the wage rate
    wage = w*eff[ij_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+r)*a[ia_com] + pen[ij_com]

    # determine labor
    if (ij_com < JR)
        lab_com = min( max( nu + (1.0-nu)*(a_plus - available)/wage, 0.0), 1.0-1e-10 )
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max( (available + wage*lab_com - a_plus) , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, ip_com, is_com] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com]

    # calculate first order condition for consumption
    return  margu(cons_com, lab_com)^(-gamma) - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, lab)

    return nu*(cons^nu*(1.0-lab)^(1.0-nu))^egam/cons

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0), 1.0-1e-10)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, is] + (1.0-varphi)*EV[ij+1, iar, ip, is], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^egam/egam + beta*psi[ij+1]*valuefunc

end 

# initializes all remaining variables
function initialize()
    global r 
    global w 
    global psi 
    global eff 
    global pen 
    global dist_theta
    global theta 
    global pi
    global eta 

    # net prices (after taxes and tranfers)
    r = 0.04
    w = 1.0

    # set survival probabilities
    psi = [1.00000, 0.99923, 0.99914, 0.99914, 0.99912, 
            0.99906, 0.99908, 0.99906, 0.99907, 0.99901,
            0.99899, 0.99896, 0.99893, 0.99890, 0.99887,
            0.99886, 0.99878, 0.99871, 0.99862, 0.99853,
            0.99841, 0.99835, 0.99819, 0.99801, 0.99785,
            0.99757, 0.99735, 0.99701, 0.99676, 0.99650,
            0.99614, 0.99581, 0.99555, 0.99503, 0.99471,
            0.99435, 0.99393, 0.99343, 0.99294, 0.99237,
            0.99190, 0.99137, 0.99085, 0.99000, 0.98871,
            0.98871, 0.98721, 0.98612, 0.98462, 0.98376,
            0.98226, 0.98062, 0.97908, 0.97682, 0.97514,
            0.97250, 0.96925, 0.96710, 0.96330, 0.95965,
            0.95619, 0.95115, 0.94677, 0.93987, 0.93445,
            0.92717, 0.91872, 0.91006, 0.90036, 0.88744,
            0.87539, 0.85936, 0.84996, 0.82889, 0.81469,
            0.79705, 0.78081, 0.76174, 0.74195, 0.72155,
            0.00000]

    # initialize age earnings process
    eff[1:JR-1] = [1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527,
            1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217,
            1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810,
            1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658,
            1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760,
            1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392,
            1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354,
            1.7701, 1.7048]

    eff[JR:JJ] .= 0.0

    # old-age transfers
    pen .= 0.0
    pen[JR:JJ] .= 0.5*sum(eff)/(JR-1)*0.33

    # initialize fixed effect
    dist_theta .= 1.0/NP
    theta[1]   = exp(-sqrt(sigma_theta))
    theta[2]   = exp( sqrt(sigma_theta))

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)

    # initialize asset grid
    grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)

    return  
end



# determines the solution to the household optimization problem
function solve_household()
    global cons_com

    # get decision in the last period of life
    for ia in 0:NA
        aplus[JJ, ia, :, :] .= 0.0
        c[JJ, ia, :, :] .= (1.0+r)*a[ia] + pen[JJ]
        l[JJ, ia, :, :] .= 0.0
        V[JJ, ia, :, :] .= valuefunc(0.0, c[JJ, ia, 1, 1], l[JJ, ia, 1, 1], JJ, 1, 1)
    end

    # interpolate individual RHS
    interpolate(JJ)

    for ij in JJ-1:-1:1

        # check about how many is to iterate
        if (ij >= JR) 
            ip_max = 1
            is_max = 1
        else
            ip_max = NP
            is_max = NS
        end

        for ia in 0:NA

            # determine decision for zero assets at retirement without pension
            if (ij >= JR && ia == 0 && pen[ij] <= 1e-10)
                aplus[ij, ia, :, :] .= 0.0
                c[ij, ia, :, :] .= 0.0
                l[ij, ia, :, :] .= 0.0
                V[ij, ia, :, :] .= valuefunc(0.0, 0.0, 0.0, ij, 1, 1)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij+1, ia, ip, is]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is

                    # solve the household problem using rootfinding
                    x_root = fzero(foc, x_in)

                    # write screen output in case of a problem
                    # check for borrowing constraint
                    if(x_root < 0.0)
                        x_root = 0.0
                        wage = w*eff[ij]*theta[ip]*eta[is]
                        available = (1.0+r)*a[ia] + pen[ij]
                        if (ij < JR)
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        cons_com = max( (available + wage*lab_com) , 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, ip, is] = x_root
                    c[ij, ia, ip, is] = cons_com
                    l[ij, ia, ip, is] = lab_com
                    V[ij, ia, ip, is] = valuefunc(x_root, cons_com, lab_com, ij, ip, is)
                end
            end

            # copy decision in retirement age
            if (ij >= JR)
                aplus[ij, ia, :, :] .= aplus[ij, ia, 1, 1]
                c[ij, ia, :, :] .= c[ij, ia, 1, 1]
                l[ij, ia, :, :] .= l[ij, ia, 1, 1]
                V[ij, ia, :, :] .= V[ij, ia, 1, 1]
            end
        end

        # interpolate individual RHS
        interpolate(ij)

        
    end
end 

# for calculating the rhs of the first order condition at age ij
function interpolate(ij)

    for ia in 0:NA
        for ip in 1:NP
            for is in 1:NS

                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is] = 0.0
                EV[ij, ia, ip, is] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ia, ip, is_p], 1e-10)
                    lhelp = max(l[ij, ia, ip, is_p], 1e-10)
                    RHS[ij, ia, ip, is] = RHS[ij, ia, ip, is] + pi[is, is_p]*margu(chelp, lhelp)
                    EV[ij, ia, ip, is]  = EV[ij, ia, ip, is] + pi[is, is_p]*V[ij, ia, ip, is_p]
                end
                RHS[ij, ia, ip, is] = ((1.0+r)*beta*psi[ij]*RHS[ij, ia, ip, is])^(-gamma)
                EV[ij, ia, ip, is] = (egam*EV[ij, ia, ip, is])^(1.0/egam)
            end
        end
    end

end


# determines the invariant distribution of households
function get_distribution()

    # set distribution to zero
    phi[:, :, :, :] .= 0.0

    # get initial distribution in age 1
    for ip in 1:NP
        phi[1, 0, ip, is_initial] = dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi = min(varphi, 1.0)

                    # redistribute households
                    for is_p in 1:NS
                        phi[ij, ial, ip, is_p] = phi[ij, ial, ip, is_p] + pi[is, is_p]*varphi*phi[ij-1, ia, ip, is]
                        phi[ij, iar, ip, is_p] = phi[ij, iar, ip, is_p] + pi[is, is_p]*(1.0-varphi)*phi[ij-1, ia, ip, is]
                    end
                end
            end
        end
    end

end 



# subroutine for calculating quantities
function aggregation()

    # calculate cohort averages
    c_coh[:] .= 0.0
    y_coh[:] .= 0.0
    l_coh[:] .= 0.0
    h_coh[:] .= 0.0
    a_coh[:] .= 0.0
    v_coh[:] .= 0.0
    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    c_coh[ij] = c_coh[ij] + c[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    y_coh[ij] = y_coh[ij] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    l_coh[ij] = l_coh[ij] + l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    h_coh[ij] = h_coh[ij] + eff[ij]*theta[ip]*eta[is]*phi[ij, ia, ip, is]
                    a_coh[ij] = a_coh[ij] + a[ia]*phi[ij, ia, ip, is]
                    v_coh[ij] = v_coh[ij] + V[ij, ia, ip, is]*phi[ij, ia, ip, is]
                end
            end
        end
    end

    # calculate cohort specific coeffcients of variation
    cv_c[:]    .= 0.0
    cv_y[:]    .= 0.0
    cv_l[:]    .= 0.0
    cv_h[:]    .= 0.0
    corr_hl[:] .= 0.0

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    cv_c[ij] = cv_c[ij] + c[ij, ia, ip, is]^2*phi[ij, ia, ip, is]
                    cv_y[ij] = cv_y[ij] + (eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is])^2*phi[ij, ia, ip, is]
                    cv_l[ij] = cv_l[ij] + l[ij, ia, ip, is]^2*phi[ij, ia, ip, is]
                    cv_h[ij] = cv_h[ij] + (eff[ij]*theta[ip]*eta[is])^2*phi[ij, ia, ip, is]
                    corr_hl[ij] = corr_hl[ij] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                end
            end
        end
    end

    #@. corr_hl = (corr_hl - h_coh*l_coh)/max(sqrt(cv_h - h_coh^2)*sqrt(cv_l - l_coh^2), 1e-10)

    @. cv_c = sqrt(cv_c - c_coh^2)/c_coh
    @. cv_y = sqrt(cv_y - y_coh^2)/max(y_coh, 1e-10)
    @. cv_l = sqrt(cv_l - l_coh^2)/max(l_coh, 1e-10)
    @. cv_h = sqrt(cv_h - h_coh^2)/max(h_coh, 1e-10)

end 
