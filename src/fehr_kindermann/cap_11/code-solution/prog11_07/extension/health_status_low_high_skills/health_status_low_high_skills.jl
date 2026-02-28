using DynamicProgrammingUtils
using Roots

# calculates year at which age ij agent is ij_p
function year(it, ij, ijp)

    year = it + ijp - ij

    if(it == 0 || year <= 0)
        year = 0
    end

    if(it == TT || year >= TT)
        year = TT
    end

    return year
end 

# function which computes the year in which the household lives
function year2(it, addit)


    year2 = it + addit

    if(year2 > TT)
        year2 = TT
    end 

    if(year2 < 0)
        year2 = 0
    end 

    if(it == 0)
        year2 = 0
    end
    if(it == TT)
        year2 = TT
    end

    return year2
end 

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


# computes the initial steady state of the economy
function get_SteadyState()


    # initialize remaining variables
    initialize()


    # iterate until value function converges
    for iter in 1:itermax

        # derive prices
        prices(0)

        # solve the household problem
        solve_household(1, 0)

        # calculate the distribution of households over state space
        get_distribution(0)

        # aggregate individual decisions
        aggregation(0)

        # determine the government parameters
        government(0)

        println(iter,"     ",round(digits = 2, HH[0]),"   ", round(digits = 2, 5.0*KK[0]/YY[0]*100.0), "   ", round(digits = 2, CC[0]/YY[0]*100.0), "   ", round(digits = 2, II[0]/YY[0]*100.0), "   ", round(digits = 2, ((1.0+r[0])^0.2-1.0)*100.0), "   ", round(digits = 2, w[0]), "   ", round(digits = 6, DIFF[0]/YY[0]*100.0))

        if(abs(DIFF[0]/YY[0])*100.0 < sig)
            break
        end

    end

end 

# initializes all remaining variables
function initialize()

    global psi 
    global eff 
    global pen 
    global dist_m
    global pi 
    global eta
    global pi_m 
    global varrho_m
    global omega
    global m
    global GAM
    global a
    global a_plus
    global dist_theta
    global theta
    global tax
    global tauc
    global tauw
    global taur
    global taup
    global kappa
    global gy
    global by
    global beq
    global BQ
    global KK
    global LL
    global YY
    global II
    global GG
    global BB
    
    global varrho
    global k
    global zeta
    global dist_zeta
    global hc

    #survival probabilities
    psi[1:6, 0, 0:TT] .= 1.00000000
    psi[7, 0, 0:TT] .= 0.98972953
    psi[8, 0, 0:TT] .= 0.98185396
    psi[9, 0, 0:TT] .= 0.97070373
    psi[10, 0, 0:TT] .= 0.95530594
    psi[11, 0, 0:TT] .= 0.93417914
    psi[12, 0, 0:TT] .= 0.90238714
    psi[13, 0, 0:TT] .= 0.83653436
    psi[14, 0, 0:TT] .= 0.71048182
    psi[15, 0, 0:TT] .= 0.52669353
    psi[16, 0, 0:TT] .= 0.31179803
    psi[17, 0, 0:TT] .= 0.00000000

    # health stauts affects survival probabilities
    psi[:, 1, 0:TT] = chi*psi[:, 0, 0:TT]
    psi[1,1,:] .= 1.0

    # initialize health shock
    dist_m[0] = 0.95
    dist_m[1] = 0.05

    # set bequest distribution
    omega[1] = 1.0/6.0
    omega[2] = 1.0/6.0
    omega[3] = 1.0/6.0
    omega[4] = 1.0/6.0
    omega[5] = 1.0/6.0
    omega[6] = 1.0/6.0
#        omega(7) = 1d0/9d0
#        omega(8) = 1d0/9d0
#        omega(9) = 1d0/9d0
    omega[7:16] .= 0.0

    # size of cohorts and bequest distribution in specific year
	for it in 0:TT
        m[1, 0, it] = dist_m[0]
        m[1, 1, it] = dist_m[1]
        itm = year(it, 2, 1)
	
		for ij in 2:JJ
			for im in 0:NM
				m[ij, im, it] = m[ij-1, im, itm]*psi[ij, im, it]/(1.0 + n_p)
			end
		end

		GAM_total = omega[1]
		
		for im in 0:NM
			for ij in 2:JJ
				GAM_total = GAM_total + omega[ij]*m[ij, im, it]
			end
		end

		for im in 0:NM
			for ij in 1:JJ
				GAM[ij, im, it] = omega[ij]/GAM_total/(NM+1)
			end
		end

	end

    # initialize age earnings process
    # initialize age earnings process
    eff[1] = 1.0000
    eff[2] = 1.3527
    eff[3] = 1.6952
    eff[4] = 1.8279
    eff[5] = 1.9606
    eff[6] = 1.9692
    eff[7] = 1.9692
    eff[8] = 1.9392
    eff[9] = 1.9007
    eff[JR:JJ] .= 0.0

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

    # get initial guess for savings decision
    for ij in 1:JJ
        for im in 0:NM
            for is in 1:NS
                @. aplus[ij, :, im, is, 0] = max(a[:]/2.0, a[1]/2.0)
            end
        end
    end


    # tax and transfers
    tax   .= 2
    tauc  .= 0.075
    tauw  .= 0.0
    taur  .= 0.0
    taup  .= 0.1
    kappa .= 0.5
    gy    = 0.19
    by    = 0.60/5.0

    beq[:, :, 0] .= 0.0
    BQ[0] = 0.0#

    # initial guesses for macro variables
    KK .= 1.0
    LL .= 1.0
    YY .= 1.0
    II .= (n_p+delta)*KK

    GG .= gy*YY[0]
    BB .= by*YY[0]

    pen .= 0.0
    pen[JR:JJ, 0] .= kappa[0]

end 

# subroutine for calculating prices
function prices(it)

    r[it] = Omega2*alpha*(KK[it]/LL[it])^(alpha-1.0)-delta
    w[it] = Omega2*(1.0-alpha)*(KK[it]/LL[it])^alpha
    rn[it] = r[it]*(1.0-taur[it])
    wn[it] = w[it]*(1.0-tauw[it]-taup[it])
    p[it] = 1.0 + tauc[it]

end 

# determines the solution to the household optimization problem
function solve_household(ij_in, it_in)
    global cons_com
    global lab_com

    # get decision in the last period of life
    it = year(it_in, ij_in, JJ)

    #bequest for the olds
    for im_ind in 0:NM
        beq[JJ, im_ind, it] = damp*GAM[JJ, im_ind, it]*BQ[it] + (1.0-damp)*beq[JJ, im_ind, it]#
    end 

    # get decision in the last period of life
    for ia in 0:NA
        aplus[JJ, ia, :, :, :, :, it] .= 0.0
        for ih in 1:NH
            c[JJ, ia, :, :, :, ih, it] .= (1.0+r)*a[ia] + pen[JJ, it] - hc[JJ, ih] + max(c_floor + hc[JJ, ih] - (1.0+r)*a[ia] - pen[JJ, it], 0.0)
            for im in 0:NM
                V[JJ, ia, :, im, :, ih, it] .= valuefunc(0.0, c[JJ, ia, 1, im, 1, ih, it], JJ, 1, im, 1)
            end
        end
    end
    
    # interpolate individual RHS
    interpolate(JJ, it)

    for ij in JJ-1:-1:1



        it = year(it_in, ij_in, ij)

        #bequest for the olds
        for im_ind in 0:NM
            beq[ij, im_ind, it] = damp*GAM[ij, im_ind, it]*BQ[it] + (1.0-damp)*beq[ij, im_ind, it]#
        end 


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
        interpolate(ij, it)

        println("Age :"*string(ij)*" DONE")
    end

end 


# for calculating the rhs of the first order condition at age ij
function interpolate(ij, it)
    
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
