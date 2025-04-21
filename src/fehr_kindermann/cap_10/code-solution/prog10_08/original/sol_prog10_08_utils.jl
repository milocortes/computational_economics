using DynamicProgrammingUtils
using Roots

# the first order condition
function foc(x_in)

    global ij_com
    global im_com
    global is_com
    global lab_com
    global cons_com

    # calculate tomorrows assets
    a_plus = x_in

    # calculate the wage rate
    wage = w*eff[ij_com]*eta[is_com]*varrho_m[im_com]

    # calculate available resources
    available = (1.0+r)*a[ia_com] + pen[ij_com]

    # determine labor
    if (ij_com < JR)
        lab_com = min(max((1.0-nu)*(a_plus-available)/wage + nu*(1.0-phi_l*Float64(im_com)), 1e-10), 1.0-phi_l*Float64(im_com)-1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max(available + wage*lab_com - a_plus, 1e-10)

    # calculate linear interpolation for future part of first order condition
    a_plus = max(a_plus, a_l)

    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, im_com, is_com] + (1.0-varphi)*RHS[ij_com+1, iar, im_com, is_com]

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, im_com)^(-gamma) - tomorrow

end 


# calculates marginal utility of consumption
function margu(cons, lab, im)

    return nu*(cons^nu*(1.0-lab-phi_l*Float64(im))^(1.0-nu))^egam/cons

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, im, is)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0), 1.0-phi_l*Float64(im)-1e-10)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, im, is] + (1.0-varphi)*EV[ij+1, iar, im, is], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help-phi_l*Float64(im))^(1.0-nu))^egam/egam + beta*psi[ij+1, im]*valuefunc

end 


# initializes all remaining variables
function initialize()
    global r 
    global w
    global psi 
    global eff 
    global pen 
    global dist_m
    global pi 
    global eta
    global pi_m 
    global varrho_m
    global a

    # net prices (after taxes and tranfers)
    r = 0.04
    w = 1.0

    # set survival probabilities
    psi[:, 0] = [1.00000, 0.99923, 0.99914, 0.99914, 0.99912,
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

    # health stauts affects survival probabilities
    psi[:, 1] = chi*psi[:, 0]

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
    pen[JR:JJ] .= 0.5*sum(eff)/Float64(JR-1)*0.33

    # initialize health shock
    dist_m[0] = 0.95
    dist_m[1] = 0.05

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)

    # probability to have bad health when current health is good
    probability_bad_health_from_good_health_temp = zeros(JJ)
    grid_Cons_Equi(probability_bad_health_from_good_health_temp, JJ, 0.1, 0.4)
    pi_m[:, 0, 1] = probability_bad_health_from_good_health_temp
    pi_m[:, 0, 0] = 1.0 .- pi_m[:, 0, 1]

    # probability to have bad health when current health is bad
    probability_bad_health_from_bad_health_temp = zeros(JJ)

    grid_Cons_Equi(probability_bad_health_from_bad_health_temp, JJ, 0.6, 0.9)
    pi_m[:, 1, 1] = probability_bad_health_from_bad_health_temp
    pi_m[:, 1, 0] = 1.0 .- pi_m[:, 1, 1]

    # initialize impact of shock on productivity
    varrho_m[0] = 1.0
    varrho_m[1] = exp(-varrho)

    # initialize asset grid
    grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)

end 


# determines the solution to the household optimization problem
function solve_household()
    global cons_com
    global lab_com

    # get decision in the last period of life
    for ia in 0:NA
        for im in 0:NM
            aplus[JJ, ia, im, :] .= 0.0
            c[JJ, ia, im, :] .= (1.0+r)*a[ia] + pen[JJ]
            l[JJ, ia, im, :] .= 0.0
            V[JJ, ia, im, :] .= valuefunc(0.0, c[JJ, ia, im, 1], l[JJ, ia, im, 1], JJ, im, 1)
        end
    end

    # interpolate individual RHS
    interpolate(JJ)

    for ij in JJ-1:-1:1

        # check about how many is to iterate
        if (ij >= JR)
            is_max = 1
        else
            is_max = NS
        end

        for ia in 0:NA
            for im in 0:NM

                # determine decision for zero assets at retirement without pension
                if (ij >= JR && ia == 0 && pen[ij] <= 1e-10)
                    aplus[ij, ia, im, :] .= 0.0
                    c[ij, ia, im, :] .= 0.0
                    l[ij, ia, im, :] .= 0.0
                    V[ij, ia, im, :] .= valuefunc(0.0, 0.0, 0.0, ij, im, 1)
                    continue
                end

                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij+1, ia, im, is]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global im_com = im
                    global is_com = is

                    # solve the household problem using rootfinding

                    x_root = fzero(foc, x_in)

                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, im, is
                    # check for borrowing constraint
                    if (x_root < 0.0)
                        x_root = 0.0
                        wage = w*eff[ij]*eta[is]*varrho_m[im]
                        available = (1.0+r)*a[ia] + pen[ij]
                        if (ij < JR)
                            lab_com = min(max(nu*(1.0-phi_l*Float64(im)) - (1.0-nu)*available/wage, 1e-10), 1.0-phi_l*Float64(im)-1e-10)
                        else
                            lab_com = 0.0
                        end
                        cons_com = max((available + wage*lab_com), 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, im, is] = x_root
                    c[ij, ia, im, is] = cons_com
                    l[ij, ia, im, is] = lab_com
                    V[ij, ia, im, is] = valuefunc(x_in, cons_com, lab_com, ij, im, is)
                end

                # copy decision in retirement age
                if (ij >= JR)
                    aplus[ij, ia, im, :] .= aplus[ij, ia, im, 1]
                    c[ij, ia, im, :] .= c[ij, ia, im, 1]
                    l[ij, ia, im, :] .= l[ij, ia, im, 1]
                    V[ij, ia, im, :] .= V[ij, ia, im, 1]
                end
            end
        end

        # interpolate individual RHS
        interpolate(ij)

        println("Age :"*string(ij)*" DONE")
    end

end 


# for calculating the rhs of the first order condition at age ij
function interpolate(ij)


    
    for ia in 0:NA
        for im in 0:NM
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ia, im, is] = 0.0
                EV[ij, ia, im, is] = 0.0
                for im_p in 0:NM
                    for is_p in 1:NS
                        chelp = max(c[ij, ia, im_p, is_p], 1e-10)
                        lhelp = max(l[ij, ia, im_p, is_p], 1e-10)
                        RHS[ij, ia, im, is] = RHS[ij, ia, im, is] + pi_m[ij, im, im_p]*pi[is, is_p]*margu(chelp, lhelp, im_p)
                        EV[ij, ia, im, is] = EV[ij, ia, im, is] + pi_m[ij, im, im_p]*pi[is, is_p]*V[ij, ia, im_p, is_p]
                    end
                end
                RHS[ij, ia, im, is] = ((1.0+r)*beta*psi[ij, im]*RHS[ij, ia, im, is])^(-gamma)
                EV[ij, ia, im, is] = (egam*EV[ij, ia, im, is])^(1.0/egam)
            end
        end
    end

end 


# determines the invariant distribution of households
function get_distribution()

    # set distribution to zero
    phi[:, :, :, :] .= 0.0

    # get initial distribution in age 1
    for im in 0:NM
        phi[1, 0, im, is_initial] = dist_m[im]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for im in 0:NM
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, im, is], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi = min(varphi, 1.0)

                    # redistribute households
                    for im_p = 0:NM
                        for is_p = 1:NS
                            phi[ij, ial, im_p, is_p] = phi[ij, ial, im_p, is_p] + pi_m[ij, im, im_p]*pi[is, is_p]*varphi*phi[ij-1, ia, im, is]
                            phi[ij, iar, im_p, is_p] = phi[ij, iar, im_p, is_p] + pi_m[ij, im, im_p]*pi[is, is_p]*(1.0-varphi)*phi[ij-1, ia, im, is]
                        end
                    end
                end
            end
        end
        println("Get Distribution Age :"*string(ij)*" DONE")
    end

end 



# subroutine for calculating quantities
function aggregation()

    # calculate fraction of good vs. bad health households
    for ij in 1:JJ
        for im in 0:NM
            frac_phi[ij, im] = sum(phi[ij, :, im, :])
        end
    end

    # calculate cohort averages
    c_coh[:, :] .= 0.0
    y_coh[:, :] .= 0.0
    a_coh[:, :] .= 0.0
    l_coh[:, :] .= 0.0
    v_coh[:, :] .= 0.0

    for ij in 1:JJ
        for ia in 0:NA
            for im in 0:NM
                for is in 1:NS
                    c_coh[ij, im]= c_coh[ij, im]+ c[ij, ia, im, is]*phi[ij, ia, im, is]/frac_phi[ij, im]
                    y_coh[ij, im]= y_coh[ij, im]+ eff[ij]*eta[is]*varrho_m[im]*phi[ij, ia, im, is]/frac_phi[ij, im]
                    a_coh[ij, im]= a_coh[ij, im]+ a[ia]*phi[ij, ia, im, is]/frac_phi[ij, im]
                    l_coh[ij, im]= l_coh[ij, im]+ l[ij, ia, im, is]*phi[ij, ia, im, is]/frac_phi[ij, im]
                    v_coh[ij, im]= v_coh[ij, im]+ V[ij, ia, im, is]*phi[ij, ia, im, is]/frac_phi[ij, im]
                end
            end
        end
    end

    # recover unconditional cohort averages
    for ij in 1:JJ
        for im in 0:NM
            c_coh[ij, NM+1] = c_coh[ij, NM+1] + c_coh[ij, im]*frac_phi[ij, im]
            y_coh[ij, NM+1] = y_coh[ij, NM+1] + y_coh[ij, im]*frac_phi[ij, im]
            a_coh[ij, NM+1] = a_coh[ij, NM+1] + a_coh[ij, im]*frac_phi[ij, im]
            l_coh[ij, NM+1] = l_coh[ij, NM+1] + l_coh[ij, im]*frac_phi[ij, im]
            v_coh[ij, NM+1] = v_coh[ij, NM+1] + v_coh[ij, im]*frac_phi[ij, im]
        end
    end

end 