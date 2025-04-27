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

# the first order condition
function foc(x_in)
    global ij_com
    global im_com
    global is_com
    global lab_com
    global cons_com
    global ia_com
    global it_com

    # calculate tomorrows assets
    a_plus = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)
    
    # get lsra transfer payment
    v_ind = v[ij_com, ia_com, im_com, is_com, it_com]

    # calculate the wage rate
    wage = wn[it_com]*eff[ij_com]*eta[is_com]*varrho_m[im_com]

    # calculate available resources
    available = (1.0+rn[it_com])*a[ia_com]+ beq[ij_com, im_com, it_com] + pen[ij_com, it_com] + v_ind

    # determine labor
    if (ij_com < JR)
        lab_com = min(max((1.0-nu)*(a_plus-available)/wage + nu*(1.0-phi_l*Float64(im_com)), 1e-10), 1.0-phi_l*Float64(im_com)-1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max((available + wage*lab_com - a_plus)/p[it_com], 1e-10)

    # calculate linear interpolation for future part of first order condition
    a_plus = max(a_plus, a_l)

    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, im_com, is_com, itp] + (1.0-varphi)*RHS[ij_com+1, iar, im_com, is_com, itp]

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, im_com, it_com)^(-gamma) - tomorrow

end 


# calculates marginal utility of consumption
function margu(cons, lab, im, it)

    return nu*(cons^nu*(1.0-lab-phi_l*Float64(im))^(1.0-nu))^egam/(p[it]*cons)

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, im, is, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0), 1.0-phi_l*Float64(im)-1e-10)

    # get tomorrows year
    itp = year(it, ij, ij+1)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, im, is, itp] + (1.0-varphi)*EV[ij+1, iar, im, is, itp], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help-phi_l*Float64(im))^(1.0-nu))^egam/egam + beta*psi[ij+1, im, itp]*valuefunc

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
        for im_ind in 0:NM
            aplus[JJ, ia, im_ind, :, it] .= 0.0
            c[JJ, ia, im_ind, :, it] .= ((1.0+rn[it])*a[ia] + beq[JJ, im_ind, it] .+ pen[JJ,it] .+ v[JJ, ia, im_ind, :, it])/p[it]
            l[JJ, ia, im_ind, :, it] .= 0.0
            V[JJ, ia, im_ind, :, it] .= valuefunc(0.0, c[JJ, ia, im_ind, 1, it], l[JJ, ia, im_ind, 1, it], JJ, im_ind, 1, it)
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
            is_max = 1
        else
            is_max = NS
        end

        for ia in 0:NA
            for im_ind in 0:NM

                # determine decision for zero assets at retirement without pension
                if (ij >= JR && ia == 0 && pen[ij] <= 1e-10)
                    aplus[ij, ia, im_ind, :] .= 0.0
                    c[ij, ia, im_ind, :] .= 0.0
                    l[ij, ia, im_ind, :] .= 0.0
                    V[ij, ia, im_ind, :] .= valuefunc(0.0, 0.0, 0.0, ij, im_ind, 1)
                    continue
                end

                for is in 1:is_max
                    # get initial guess for the individual choices
                    x_in = aplus[ij+1, ia, im_ind, is, it]

                    println(ij, ",", is, ",", ia, ",", im_ind, ",", x_in)

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global im_com = im_ind
                    global is_com = is
                    global it_com = it

                    # solve the household problem using rootfinding

                    x_root = find_zero(foc, x_in)

                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, im, is
                    # check for borrowing constraint
                    if (x_root < 0.0)
                        x_root = 0.0
                        wage = wn[it]*eff[ij]*eta[is]*varrho_m[im_ind]
                        v_ind = v[ij, ia, im_ind, is, it]
                        available = (1.0+rn[it])*a[ia] + beq_coh[ij, im_ind, it] + pen[ij, it]  
                        if (ij < JR)
                            global lab_com = min(max(nu*(1.0-phi_l*Float64(im_ind)) - (1.0-nu)*available/wage, 1e-10), 1.0-phi_l*Float64(im_ind)-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        global cons_com = max((available + wage*lab_com), 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, im_ind, is, it] = x_root
                    c[ij, ia, im_ind, is, it] = cons_com
                    l[ij, ia, im_ind, is, it] = lab_com
                    V[ij, ia, im_ind, is, it] = valuefunc(x_in, cons_com, lab_com, ij, im_ind, is, it)
                end

                # copy decision in retirement age
                if (ij >= JR)
                    aplus[ij, ia, im_ind, :, it] .= aplus[ij, ia, im_ind, 1, it]
                    c[ij, ia, im_ind, :, it] .= c[ij, ia, im_ind, 1, it]
                    l[ij, ia, im_ind, :, it] .= l[ij, ia, im_ind, 1, it]
                    V[ij, ia, im_ind, :, it] .= V[ij, ia, im_ind, 1, it]
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
        for im in 0:NM
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ia, im, is, it] = 0.0
                EV[ij, ia, im, is, it] = 0.0
                for im_p in 0:NM
                    for is_p in 1:NS
                        chelp = max(c[ij, ia, im_p, is_p, it], 1e-10)
                        lhelp = max(l[ij, ia, im_p, is_p, it], 1e-10)
                        RHS[ij, ia, im, is, it] = RHS[ij, ia, im, is, it] + pi_m[ij, im, im_p]*pi[is, is_p]*margu(chelp, lhelp, im_p, it)
                        EV[ij, ia, im, is, it] = EV[ij, ia, im, is, it] + pi_m[ij, im, im_p]*pi[is, is_p]*V[ij, ia, im_p, is_p, it]
                    end
                end
                RHS[ij, ia, im, is, it] = ((1.0+rn[it])*beta*psi[ij, im, it]*RHS[ij, ia, im, is, it])^(-gamma)
                EV[ij, ia, im, is, it] = (egam*EV[ij, ia, im, is, it])^(1.0/egam)
            end
        end
    end

end 