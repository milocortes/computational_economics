using DynamicProgrammingUtils
using Roots


#  the first order condition
function foc(x_in)

    global ij_com
    global ip_com
    global im_com
    global is_com
    global cons_com
    global ih_com

    #  calculate tomorrows assets
    a_plus = x_in

    #  calculate the wage rate
    wage = w*eff[ij_com]*theta[ip_com]*eta[is_com]*varrho_m[ip_com, im_com]

    #  calculate available resources
    btran = max(c_floor + hc[ij_com, ih_com] - (1.0+r)*(a[ia_com] + a_bor[ij_com, ip_com]) - wage - pen[ij_com], 0.0)
    available = (1.0+r)*(a[ia_com] + a_bor[ij_com, ip_com]) + wage + pen[ij_com] - hc[ij_com, ih_com] + btran

    #  calculate consumption
    cons_com = available - (a_plus + a_bor[ij_com+1, ip_com])

    #  calculate linear interpolation for future part of foc
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, ip_com, im_com, is_com] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, im_com, is_com]

    #  calculate first order condition for consumption
    return cons_com - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, im)
    return (1.0-delta*im)*cons^(-1.0/gamma)
end


# calculates the value function
function valuefunc(a_plus, cons, ij, ip, im, is)


    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, im, is] + (1.0-varphi)*EV[ij+1, iar, ip, im, is], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (1.0-delta*im)*c_help^egam/egam + beta*psi[ij+1, im]*valuefunc

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
    global pi_m
    global varrho
    global varrho_m
    global k
    global zeta
    global dist_zeta
    global hc

    # net prices (after taxes and tranfers)
    r = 0.01
    w = 1.0

    # set survival probabilities
    psi[:, 0] .= [1.0, 0.9884, 0.9977, 0.9977, 0.9989, 
    0.9989, 0.9989, 0.9989, 0.9989, 0.9985, 
    0.9985, 0.9985, 0.9985, 0.9985, 0.9961, 
    0.9961, 0.9961, 0.9961, 0.9961, 0.9941, 
    0.9941, 0.9941, 0.9941, 0.9941, 0.9935, 
    0.9935, 0.9935, 0.9935, 0.9935, 0.992, 
    0.992, 0.992, 0.992, 0.992, 0.9888, 
    0.9888, 0.9888, 0.9888, 0.9888, 0.9863, 
    0.9863, 0.9863, 0.9863, 0.9863, 0.9808, 
    0.9808, 0.9808, 0.9808, 0.9808, 0.9716, 
    0.9716, 0.9716, 0.9716, 0.9716, 0.9574, 
    0.9574, 0.9574, 0.9574, 0.9574, 0.9368, 
    0.9368, 0.9368, 0.9368, 0.9368, 0.9068, 
    0.9068, 0.9068, 0.9068, 0.9068, 0.8627, 
    0.8627, 0.8627, 0.8627, 0.8627, 0.7911, 
    0.7911, 0.7911, 0.7911, 0.7911, 0.7911, 
    0.0]

    # health depended survival probabilities
    psi[:, 1] .= chi*psi[:, 0]

    # initialize age earnings process
    eff[1:JR-1] = [1.0, 1.02120569292179, 1.09455534495088, 1.12860287958448, 1.22638160666327, 1.23391086554654, 
        1.30222364397719, 1.30147188600755, 1.32093252269188, 1.3350357635024, 1.34108586627201, 1.35264098907184,
        1.42796373993807, 1.39899505168983, 1.3725561241996, 1.27471012406845, 1.348005595511, 1.33813911364699, 
        1.3270604580725, 1.3992304207336, 1.29502534677326, 1.42368629438838, 1.41635229155533, 1.39785527090306,
        1.33458108162614, 1.32920971541425, 1.31784385225717, 1.24262648766814, 1.31066796168053, 1.26712374306502,
        1.26233154950005, 1.2475788827581, 1.29695927207029, 1.24823644891577, 1.26796037783097, 1.19068172151984, 
        1.30560560553093, 1.34965626605632, 1.11725915571105, 1.29358505316261, 1.15314506883214, 1.15961363869665,
        1.16021484018955, 0.985727573756957]
    eff[JR:JJ] .= 0.0

    # old-age transfers
    pen .= 0.0
    pen[JR:JJ] .= 0.62*w*sum(eff)/Float64(JR-1)

    # initialize fixed effect
    dist_theta .= 1.0/Float64(NP)
    theta[1]   = exp(-sqrt(sigma_theta))
    theta[2]   = exp( sqrt(sigma_theta))

    # initialize health shock
    dist_m[0] = 0.95
    dist_m[1] = 0.05

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)


    # probability to have bad health when current health is good
    prob_bad_health_from_good_1 = zeros(JJ)
    prob_bad_health_from_good_2 = zeros(JJ)

    grid_Cons_Equi(prob_bad_health_from_good_1, JJ, 0.1, 0.4)
    grid_Cons_Equi(prob_bad_health_from_good_2, JJ, 0.0, 0.2)

    pi_m[:, 1, 0, 1] = prob_bad_health_from_good_1
    pi_m[:, 2, 0, 1] = prob_bad_health_from_good_2

    pi_m[:, :, 0, 0] = 1.0 .- pi_m[:, :, 0, 1]

    # probability to have bad health when current health is bad
    prob_bad_health_from_bad_1 = zeros(JJ)
    prob_bad_health_from_bad_2 = zeros(JJ)

    grid_Cons_Equi(prob_bad_health_from_bad_1, JJ, 0.6, 0.9)
    grid_Cons_Equi(prob_bad_health_from_bad_2, JJ, 0.6, 0.9)

    pi_m[:, 1, 1, 1] = prob_bad_health_from_bad_1
    pi_m[:, 2, 1, 1] = prob_bad_health_from_bad_2

    pi_m[:, :, 1, 0] = 1.0 .- pi_m[:, :, 1, 1]

    # initialize impact of shock on productivity
    varrho[1] = 0.2
    varrho[2] = 0.1
    for ip in 1:NP
        for im in 0:NM
            varrho_m[ip, im] = exp(-varrho[ip]*Float64(im))
        end
    end

    # calculate the out-of-pocket expenses for health
    grid_Cons_Equi(k, JJ, 0.3, 0.9)

    # normally distributed stochastic term
    normal_discrete(zeta, dist_zeta, NH, 0.0, sigma_zeta)
    zeta = exp.(zeta)

    # out of pocket expenses
    for ij in 1:JJ
        hc[ij, :] = k[ij]*zeta
    end

    # initialize asset grid
    grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)

    # value of the borrowing constraint
    a_bor .= 0.0

    # calculate endogenous borrowing constraints
    for ij in 1:JJ
        for ip in 1:NP

            # calculate natural borrowing limit
            abor_temp = 0.0
            for ijj in JJ:-ij:1
                
                abor_temp = abor_temp/(1.0+r) + eff[ijj]*theta[ip]*eta[1] + pen[ijj]
            end
            
            abor_temp = min(-abor_temp/(1.0+r)+1e-4, 0.0)
            # set maximum of natural and exogenous borrowing constraint
            a_bor[ij, ip] = max(a_bor[ij, ip], abor_temp)
        end
    end

end 


# determines the solution to the household optimization problem
function solve_household()
    global cons_com

    # get decision in the last period of life
    for ia in 0:NA
        aplus[JJ, ia, :, :, :, :] .= 0.0
        for ih in 1:NH
            c[JJ, ia, :, :, :, ih] .= (1.0+r)*a[ia] + pen[JJ] - hc[JJ, ih] + max(c_floor + hc[JJ, ih] - (1.0+r)*a[ia] - pen[JJ], 0.0)
            for im in 0:NM
                V[JJ, ia, :, im, :, ih] .= valuefunc(0.0, c[JJ, ia, 1, im, 1, ih], JJ, 1, im, 1)
            end
        end
    end

    # interpolate individual RHS
    interpolate(JJ)

    for ij = JJ-1:-1:1

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
                aplus[ij, ia, :, :, :, :] .= 0.0
                c[ij, ia, :, :, :, :] .= 0.0
                for im in 0:NM
                    V[ij, ia, :, im, :, :] .= valuefunc(0.0, 0.0, ij, 1, im, 1)
                end
                continue
            end

            # solve household problem using rootfinding method
            for ip in 1:ip_max
                for im in 0:NM
                    for is in 1:is_max
                        for ih in 1:NH

                            # get initial guess for the individual choices
                            x_in = aplus[ij+1, ia, ip, im, is, ih]

                            # set up communication variables
                            global ij_com = ij
                            global ia_com = ia
                            global ip_com = ip
                            global im_com = im
                            global is_com = is
                            global ih_com = ih

                            # solve the household problem using rootfinding
                            x_root = fzero(foc, x_in)

                            # write screen output in case of a problem
                            #if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, im, is, ih

                            # check for borrowing constraint
                            if (x_root < 0.0)
                                x_root = 0.0
                                wage = w*eff[ij]*theta[ip]*eta[is]*varrho_m[ip, im]
                                global cons_com = max((1.0+r)*a[ia] + wage + pen[ij] - hc[ij, ih] + max(c_floor + hc[ij, ih] - (1.0+r)*a[ia] - wage - pen[ij], 0.0), 1e-10)
                            end

                            # copy decisions
                            aplus[ij, ia, ip, im, is, ih] = x_root
                            c[ij, ia, ip, im, is, ih] = cons_com
                            V[ij, ia, ip, im, is, ih] = valuefunc(x_root, cons_com, ij, ip, im, is)
                        end
                    end
                end
            end

            # copy decision in retirement age
            if (ij >= JR)
                for im in 0:NM
                    for ih in 1:NH
                        aplus[ij, ia, :, im, :, ih] .= aplus[ij, ia, 1, im, 1, ih]
                        c[ij, ia, :, im, :, ih] .= c[ij, ia, 1, im, 1, ih]
                        V[ij, ia, :, im, :, ih] .= V[ij, ia, 1, im, 1, ih]
                    end
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
        for ip in 1:NP
            for im in 0:NM
                for is in 1:NS

                    # calculate the RHS of the first order condition
                    RHS[ij, ia, ip, im, is] = 0.0
                    EV[ij, ia, ip, im, is] = 0.0
                    for im_p in 0:NM
                        for is_p in 1:NS
                            for ih_p in 1:NH
                                chelp = max(c[ij, ia, ip, im_p, is_p, ih_p], 1e-10)
                                RHS[ij, ia, ip, im, is] = RHS[ij, ia, ip, im, is] + dist_zeta[ih_p]* pi_m[ij, ip, im, im_p]*pi[is, is_p]*margu(chelp, im_p)
                                EV[ij, ia, ip, im, is] = EV[ij, ia, ip, im, is] + dist_zeta[ih_p]*pi_m[ij, ip, im, im_p]*pi[is, is_p]*V[ij, ia, ip, im_p, is_p, ih_p]
                            end
                        end
                    end
                    RHS[ij, ia, ip, im, is] = (beta*psi[ij, im]*(1.0+r)/(1.0-delta*Float64(im))*RHS[ij, ia, ip, im, is])^(-gamma)
                    EV[ij, ia, ip, im, is] = (egam*EV[ij, ia, ip, im, is])^(1.0/egam)
                end
            end
        end
    end

    return
end 


# determines the invariant distribution of households
function get_distribution()

    # set distribution to zero
    phi[:, :, :, :, :, :] .= 0.0;

    # get initial distribution in age 1
    for ip in 1:NP
        for im in 0:NM
            for ih in 1:NH
                
                # find the zero on the asset grid
                ial, iar, varphi = linint_Grow(-a_bor[1, ip], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                phi[1, ial, ip, im, is_initial, ih] = varphi*dist_theta[ip]*dist_m[im]*dist_zeta[ih]
                phi[1, iar, ip, im, is_initial, ih] = (1.0-varphi)*dist_theta[ip]*dist_m[im]*dist_zeta[ih]

            end
        end
    end

    # successively compute distribution over ages
    for ij in 2:JJ
        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for im in 0:NM
                    for is in 1:NS
                        for ih in 1:NH

                            # interpolate yesterday's savings decision
                            ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, im, is, ih], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                            # restrict values to grid just in case
                            ial = min(ial, NA)
                            iar = min(iar, NA)
                            varphi = min(varphi, 1.0)

                            # redistribute households
                            for im_p in 0:NM
                                for is_p in 1:NS
                                    for ih_p in 1:NH
                                        phi[ij, ial, ip, im_p, is_p, ih_p] = phi[ij, ial, ip, im_p, is_p, ih_p] + pi[is, is_p]*pi_m[ij, ip, im, im_p]*varphi*dist_zeta[ih_p]*phi[ij-1, ia, ip, im, is, ih]
                                        phi[ij, iar, ip, im_p, is_p, ih_p] = phi[ij, iar, ip, im_p, is_p, ih_p] + pi[is, is_p]*pi_m[ij, ip, im, im_p]*(1.0-varphi)*dist_zeta[ih_p]*phi[ij-1, ia, ip, im, is, ih]
                                    end
                                end
                            end
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
        for ip in 1:NP
            for im in 0:NM
                frac_phi[ij, ip, im] = sum(phi[ij, :, ip, im, :, :])
            end
        end
    end

    # calculate cohort averages
    c_coh[:, :, :] .= 0.0
    y_coh[:, :, :] .= 0.0
    a_coh[:, :, :] .= 0.0
    v_coh[:, :, :] .= 0.0
    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for im in 0:NM
                    for is in 1:NS
                        for ih in 1:NH
                            c_coh[ij, ip, im] = c_coh[ij, ip, im] + c[ij, ia, ip, im, is, ih]*phi[ij, ia, ip, im, is, ih]/frac_phi[ij, ip, im]
                            y_coh[ij, ip, im] = y_coh[ij, ip, im]+ eff[ij]*theta[ip]*eta[is]^varrho_m[ip, im]*phi[ij, ia, ip, im, is, ih]/frac_phi[ij, ip, im]
                            a_coh[ij, ip, im] = a_coh[ij, ip, im]+ a[ia]*phi[ij, ia, ip, im, is, ih]/frac_phi[ij, ip, im]
                            v_coh[ij, ip, im] = v_coh[ij, ip, im]+ V[ij, ia, ip, im, is, ih]*phi[ij, ia, ip, im, is, ih]/frac_phi[ij, ip, im]
                        end
                    end
                end
            end
        end
    end

    # recover unconditional cohort averages
    for ij = 1:JJ
        for ip = 1:NP
            for im = 0:NM
                c_coh[ij, ip, NM+1] = c_coh[ij, ip, NM+1] + c_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, ip, :])
                y_coh[ij, ip, NM+1] = y_coh[ij, ip, NM+1] + y_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, ip, :])
                a_coh[ij, ip, NM+1] = a_coh[ij, ip, NM+1] + a_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, ip, :])
                v_coh[ij, ip, NM+1] = v_coh[ij, ip, NM+1] + v_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, ip, :])
                c_coh[ij, NP+1, im] = c_coh[ij, NP+1, im] + c_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, im])
                y_coh[ij, NP+1, im] = y_coh[ij, NP+1, im] + y_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, im])
                a_coh[ij, NP+1, im] = a_coh[ij, NP+1, im] + a_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, im])
                v_coh[ij, NP+1, im] = v_coh[ij, NP+1, im] + v_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, im])
                c_coh[ij, NP+1, NM+1] = c_coh[ij, NP+1, NM+1] + c_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, :])
                y_coh[ij, NP+1, NM+1] = y_coh[ij, NP+1, NM+1] + y_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, :])
                a_coh[ij, NP+1, NM+1] = a_coh[ij, NP+1, NM+1] + a_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, :])
                v_coh[ij, NP+1, NM+1] = v_coh[ij, NP+1, NM+1] + v_coh[ij, ip, im]*frac_phi[ij, ip, im]/sum(frac_phi[ij, :, :])
            end
        end
    end

end 
