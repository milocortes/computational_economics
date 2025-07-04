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
    global ia_com
    global ip_com
    global is_com
    global it_com
    global cons_com
    global lab_com
    global ik_com

    # calculate tomorrows assets
    a_plus  = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)

    # get lsra transfer payment
    v_ind = v[ij_com, ik_com, ia_com, ip_com, is_com, it_com]

    # calculate the wage rate
    wage = wn[it_com]*eff[ij_com,ik_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    #available = (1.0+rn[it_com])*a[ia_com] + pen[ij_com, it_com] + v_ind
    available = (1.0+rn[it_com])*a[ia_com] + beq[ij_com, ik_com, it_com] + pen[ij_com, it_com] + v_ind

    # determine labor
    if (ij_com < JR)
        lab_com = min( max( nu + (1.0-nu)*(a_plus-available)/wage, 0.0) , 1.0-1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max( (available + wage*lab_com - a_plus)/p[it_com] , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = max(varphi*RHS[ij_com+1, ik_com, ial, ip_com, is_com, itp] + (1.0-varphi)*RHS[ij_com+1, ik_com, iar, ip_com, is_com, itp], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, it_com)^(-gamma) - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, lab, it)

    margu = nu*(cons^nu*(1.0-lab)^(1.0-nu))^(1.0-1.0/gamma)/(p[it]*cons)

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ik, ip, is, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0),1.0-1e-10)

    # get tomorrows year
    itp = year(it, ij, ij+1)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0

    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ik, ial, ip, is, itp] + (1.0-varphi)*EV[ij+1, ik, iar, ip, is, itp], 1e-10)^(1.0-1.0/gamma)/(1.0-1.0/gamma)
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^(1.0-1.0/gamma)/(1.0-1.0/gamma) + beta*psi[ij+1, itp]valuefunc

end 



# computes the initial steady state of the economy
function get_SteadyState()

    # initialize remaining variables
    initialize()

    # start timer
    #call tic()

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

    #call toc
    #call output(0)

    #write(*,'(a/)')'ATTENTION: NO CONVERGENCE ###'

end



# initializes the remaining model parameters and variables
function initialize()
    global psi
    global omega
    global m
    global GAM
    global a
    global a_plus
    global eff
    global dist_theta
    global theta
    global pi 
    global eta
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
    global pen
    global pop 

    println("INITIAL EQUILIBRIUM")
    println("ITER     H     K/Y     C/Y     I/Y       r       w        DIFF")

    #= set up population structure
    for ij in 1:JJ
        pop[ij, 0] = 1.0/(1.0+n_p)^(ij-1)
    end
    for ij in 1:JJ
        m[ij, 0] = pop[ij, 0]/pop[1, 0]
    end

    for it in 0:TT
        m[1, 1, it] = 1 - 0.6
        m[1, 2, it] = 1.0 - m[1, 1, it]
        
        itm = year2(it, -1)

        for ik in 1:SS
            for ij in 2:JJ
                m[ij, ik, it] = m[ij-1, ik, itm]/(1.0+n_p)
            end
        end
    end
    =#

    #survival probabilities
    psi[1:6,0:TT] .= 1.00000000
    psi[7,0:TT] .= 0.98972953
    psi[8,0:TT] .= 0.98185396
    psi[9,0:TT] .= 0.97070373
    psi[10,0:TT] .= 0.95530594
    psi[11,0:TT] .= 0.93417914
    psi[12,0:TT] .= 0.90238714
    psi[13,0:TT] .= 0.83653436
    psi[14,0:TT] .= 0.71048182
    psi[15,0:TT] .= 0.52669353
    psi[16,0:TT] .= 0.31179803
    psi[17,0:TT] .= 0.00000000

#        psi(:,0:TT) = 1.0d0

    # set bequest distribution
    omega[1] = 1.0/9.0
    omega[2] = 1.0/9.0
    omega[3] = 1.0/9.0
    omega[4] = 1.0/9.0
    omega[5] = 1.0/9.0
    omega[6] = 1.0/9.0
    omega[7] = 1.0/9.0
    omega[8] = 1.0/9.0
    omega[9] = 1.0/9.0
    omega[10:16] .= 0.0

	# size of cohorts and bequest distribution in specific year
	for it in 0:TT
        m[1, 1, it] = 1.0 - 0.6
        m[1, 2, it] = 1.0 - m[1, 1, it]

        itm = year(it, 2, 1)
	
		for ij in 2:JJ
			for ik in 1:SS
				m[ij, ik, it] = m[ij-1, ik, itm]*psi[ij, it]/(1.0 + n_p)
			end
		end

		GAM_total = omega[1]
		
		for ik in 1:SS
			for ij in 2:JJ
				GAM_total = GAM_total + omega[ij]*m[ij, ik, it]
			end
		end

		for ik in 1:SS
			for ij in 1:JJ
				GAM[ij, ik, it] = omega[ij]/GAM_total/SS
			end
		end
	end

    # initialize asset grid
    grid_Cons_Grow(a,  NA+1, a_l, a_u, a_grow)

    # get initial guess for savings decision
    for ij in 1:JJ
        for ik in 1:SS
            for ip in 1:NP
                for is in 1:NS
                    @. aplus[ij,ik, :, ip, is, 0] = max(a[:]/2.0, a[1]/2.0)
                end
            end
        end
    end

    # initialize age earnings process
    eff[1,2] = 1.0000
    eff[2,2] = 1.3527
    eff[3,2] = 1.6952
    eff[4,2] = 1.8279
    eff[5,2] = 1.9606
    eff[6,2] = 1.9692
    eff[7,2] = 1.9692
    eff[8,2] = 1.9392
    eff[9,2] = 1.9007
    eff[JR:JJ,2] .= 0.0

    eff[:,1] = eff[:,2]*0.8

    # initialize fixed effect
    dist_theta .= 1.0/Float64(NP)
    theta[1]   = -sqrt(sigma_theta)
    theta[2]   = sqrt(sigma_theta)
    theta = exp.(theta)

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)

    # tax and transfers
    tax   .= 2
    tauc  .= 0.075
    tauw  .= 0.0
    taur  .= 0.0
    taup  .= 0.1
    kappa .= 0.5
    gy    = 0.19
    by    = 0.60/5.0

    beq[:,:, 0] .= 0.0
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



# function for calculating prices
function prices(it)

    r[it] = Omega*alpha*(KK[it]/LL[it])^(alpha-1.0)-delta
    w[it] = Omega*(1.0-alpha)*(KK[it]/LL[it])^alpha
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
    for ik in 1:SS
       beq[JJ, ik, it] = damp*GAM[JJ, ik, it]*BQ[it] + (1.0-damp)*beq[JJ, ik, it]#
    end

    for ik in 1:SS
        for ia in 0:NA
            aplus[JJ, ik, ia, :, :, it] .= 0.0
            #c[JJ, ik, ia, :, :, it] .= ((1.0+rn[it])*a[ia] .+ pen[JJ, it] .+ v[JJ, ik, ia, :, :, it])/p[it]
            c[JJ, ik, ia, :, :, it] .= ((1.0+rn[it])*a[ia]+ beq[JJ, ik, it] .+ pen[JJ, it] .+ v[JJ, ik, ia, :, :, it])/p[it]#
            l[JJ, ik, ia, :, :, it] .= 0.0
            VV[JJ, ik, ia, :, :, it] .= valuefunc(0.0, c[JJ, ik, ia, 1, 1, it],l[JJ, ik, ia, 1, 1, it], JJ, ik, 1, 1, it)
        end
    end

    # interpolate individual RHS
    interpolate(JJ, 1, it)
    interpolate(JJ, 2, it)

    for ij in JJ-1:-1:ij_in
        
        it = year(it_in, ij_in, ij)

        for ik in 1:SS

        #bequest for the olds
        beq[ij, ik, it] = damp*GAM[ij, ik, it]*BQ[it] + (1.0-damp)*beq[ij, ik, it]#

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
            if (ij >= JR && ia == 0 && kappa[it] <= 1e-10)
                aplus[ij, ik, ia, :, :, it] .= 0.0
                c[ij, ik, ia, :, :, it] .= 0.0
                l[ij, ik, ia, :, :, it] .= 0.0
                VV[ij, ik, ia, :, :, it] .= valuefunc(0.0, 0.0, 0.0, ij, ik, 1, 1, it)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij, ik, ia, ip, is, it]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is
                    global it_com = it
                    global ik_com = ik

                    # solve the household problem using rootfinding
                    #call fzero(x_in, foc, check)
                    x_root = fzero(foc, x_in)

                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is, it

                    # check for borrowing constraint
                    if (x_root < 0.0)
                        x_root = 0.0
                        wage = wn[it]*eff[ij, ik]*theta[ip]*eta[is]
                        v_ind = v[ij, ik, ia, ip, is, it]
                        available = (1.0+rn[it])*a[ia] + beq[ij, ik, it] + pen[ij, it] + v_ind
                        if (ij < JR)
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        global cons_com = max( (available + wage*lab_com)/p[it] , 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ik, ia, ip, is, it] = x_root
                    c[ij, ik, ia, ip, is, it] = cons_com
                    l[ij, ik, ia, ip, is, it] = lab_com
                    VV[ij, ik, ia, ip, is, it] = valuefunc(x_root, cons_com, lab_com, ij, ik, ip, is, it)

                end

                # copy decision in retirement age
                if (ij >= JR)
                    aplus[ij, ik, ia, :, :, it] .= aplus[ij, ik, ia, 1, 1, it]
                    c[ij, ik, ia, :, :, it] .= c[ij, ik, ia, 1, 1, it]
                    l[ij, ik, ia, :, :, it] .= l[ij, ik, ia, 1, 1, it]
                    VV[ij, ik, ia, :, :, it] .= VV[ij, ik, ia, 1, 1, it]
                end
            end
        end
        # interpolate individual RHS
        interpolate(ij, ik, it)

        end
    end

end


# for calculating the rhs of the first order condition at age ij
function interpolate(ij, ik, it)

    for ia in 0:NA
        for ip in 1:NP
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ik, ia, ip, is, it] = 0.0
                EV[ij, ik, ia, ip, is, it] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ik, ia, ip, is_p, it],1e-10)
                    lhelp = max(l[ij, ik, ia, ip, is_p, it],1e-10)
                    RHS[ij, ik, ia, ip, is, it] = RHS[ij, ik, ia, ip, is, it] + pi[is, is_p]*margu(chelp, lhelp, it)
                    EV[ij, ik, ia, ip, is, it]  = EV[ij, ik, ia, ip, is, it] + pi[is, is_p]*VV[ij, ik, ia, ip, is_p, it]
                end
                #RHS[ij, ik, ia, ip, is, it] = ((1.0+rn[it])*beta*RHS[ij, ik, ia, ip, is, it])^(-gamma)
                #RHS[ij, ik, ia, ip, is, it] = ((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ik, ia, ip, is, it])^(-gamma)
                RHS[ij, ik, ia, ip, is, it] = max((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ik, ia, ip, is, it],1e-10)^(-gamma)
                EV[ij, ik, ia, ip, is, it] = ((1.0-1.0/gamma)*EV[ij, ik, ia, ip, is, it])^(1.0/(1.0-1.0/gamma))
            end
        end
    end

end


# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, :, it] .= 0.0

    # get initial distribution in age 
    for ik in 1:SS
        for ip in 1:NP
            phi[1, ik, 0, ip, is_initial, it] = dist_theta[ip]
        end
    end

    # successively compute distribution over ages
    for ij in 2:JJ
        for ik in 1:SS
            # iterate over yesterdays gridpoints
            for ia in 0:NA
                for ip in 1:NP
                    for is in 1:NS

                        # interpolate yesterday's savings decision
                        ial, iar, varphi = linint_Grow(aplus[ij-1, ik, ia, ip, is, itm], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                        # restrict values to grid just in case
                        ial = max(min(ial, NA-1), 0)
                        iar = max(min(iar, NA), 1)
                        varphi = max(min(varphi, 1.0), 0.0)

                        # redistribute households
                        for is_p in 1:NS
                            phi[ij, ik, ial, ip, is_p, it] = phi[ij, ik, ial, ip, is_p, it] + pi[is, is_p]*varphi*phi[ij-1, ik, ia, ip, is, itm]
                            phi[ij, ik, iar, ip, is_p, it] = phi[ij, ik, iar, ip, is_p, it] + pi[is,is_p]*(1.0-varphi)*phi[ij-1, ik, ia,ip,is,itm]
                        end
                    end
                end
            end
        end
    end

end


# function for calculating quantities in a certain
function aggregation(it)

    # get tomorrow's year
    itp = year(it, 1, 2)
    LL_old = LL[it]

    m_coh = zeros(JJ, SS)

    # calculate cohort aggregates
    c_coh[:, :, it]  .= 0.0
    l_coh[:, :, it]  .= 0.0
    y_coh[:, :, it]  .= 0.0
    a_coh[:, :, it]  .= 0.0
    VV_coh[:, :, it] .= 0.0
    m_coh[:, :]      .= 0.0
    FLC[:, :, it]     .= 0.0
    beq_coh[:,:, it] .= 0.0#

    for ij in 1:JJ
        for ik in 1:SS
            for ia in 0:NA
                for ip in 1:NP
                    for is in 1:NS
                        c_coh[ij, ik, it] = c_coh[ij, ik, it] + c[ij, ik, ia, ip, is, it]*phi[ij, ik, ia, ip, is, it]
                        l_coh[ij, ik, it] = l_coh[ij, ik, it] + l[ij, ik, ia, ip, is, it]*phi[ij, ik, ia, ip, is, it]
                        y_coh[ij, ik, it] = y_coh[ij, ik, it] + eff[ij, ik]*theta[ip]*eta[is]*l[ij, ik, ia, ip, is, it]*phi[ij, ik, ia, ip, is, it]
                        a_coh[ij, ik, it] = a_coh[ij, ik, it] + a[ia]*phi[ij, ik, ia, ip, is, it]

                        # exclude households who die
                        if(ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                            continue
                        end
                        if(aplus[ij, ik, ia, ip, is, it] < 1e-4)
                            FLC[ij, ik, it] = FLC[ij, ik, it] + phi[ij, ik, ia, ip, is, it]
                        end
                        VV_coh[ij, ik, it] = VV_coh[ij, ik, it] + VV[ij, ik, ia, ip, is, it]*phi[ij, ik, ia, ip, is, it]
                        m_coh[ij, ik]      = m_coh[ij, ik] + phi[ij, ik, ia, ip, is, it]
                        beq_coh[ij, ik, it] = beq_coh[ij, ik, it] + a[ia]*(1.0 + rn[it])*(1.0 - psi[ij, it])*phi[ij, ik, ia, ip, is, it]/psi[ij, it] #

                    end
                end
            end
        end
    end

    # normalize VV_coh (because hh excluded)
    VV_coh[:, :, it] = VV_coh[:, :, it]./m_coh
    FLC[:, :, it] = FLC[:, :, it]./m_coh

    # calculate aggregate quantities
    CC[it] = 0.0
    LL[it] = 0.0
    HH[it] = 0.0
    AA[it] = 0.0
    workpop = 0.0
    BQ[it] = 0.0

    for ij in 1:JJ
        for ik in 1:SS
            CC[it] = CC[it] + c_coh[ij, ik, it]*m[ij, ik, it]
            LL[it] = LL[it] + y_coh[ij, ik, it]*m[ij, ik, it]
            HH[it] = HH[it] + l_coh[ij, ik, it]*m[ij, ik, it]
            #AA[it] = AA[it] + a_coh[ij, ik, it]*m[ij, ik, it]
            AA[it] = AA[it] + a_coh[ij, ik, it]*m[ij, ik, it]/psi[ij, it]#
            BQ[it] = BQ[it] + beq_coh[ij, ik, it]*m[ij, ik, it]   #
            if (ij < JR)
                workpop = workpop + m[ij, ik, it]
            end
        end
    end

    # damping and other quantities
    KK[it] = damp*(AA[it]-BB[it]-BA[it])+(1.0-damp)*KK[it]
    LL[it] = damp*LL[it] + (1.0-damp)*LL_old
    II[it] = (1.0+n_p)*KK[itp] - (1.0-delta)*KK[it]
    YY[it] = Omega * KK[it]^alpha * LL[it]^(1.0-alpha)

    # get average income and average working hours
    INC[it] = w[it]*LL[it]/workpop
    HH[it]  = HH[it]/workpop

    # get difference on goods market
    DIFF[it] = YY[it]-CC[it]-II[it]-GG[it]

end


# function for calculating government parameters
function government(it)

    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # set government quantities and pension payments
    GG[it] = gy*YY[0]
    BB[it] = by*YY[0]
    pen[JR:JJ, it] .= kappa[it]*INC[itm]
    PP[it] = 0.0

    for ij in JR:JJ
        for ik in 1:SS
            PP[it] = PP[it] + pen[ij, it]*m[ij, ik, it]
        end
    end

    # calculate government expenditure
    expend = GG[it] + (1.0+r[it])*BB[it] - (1.0+n_p)*BB[itp]

    # get budget balancing tax rate
    if (tax[it] == 1)
        tauc[it] = (expend - (tauw[it]*w[it]*LL[it] + taur[it]*r[it]*AA[it]))/CC[it]
        p[it]    = 1.0 + tauc[it]
    elseif (tax[it] == 2)
        tauw[it] = (expend - tauc[it]*CC[it])/(w[it]*LL[it] + r[it]*AA[it])
        taur[it] = tauw[it]
    elseif (tax[it] == 3)
        tauw[it] = (expend - (tauc[it]*CC[it] + taur[it]*r[it]*AA[it]))/(w[it]*LL[it])
    else
        taur[it] = (expend - (tauc[it]*CC[it] + tauw[it]*w[it]*LL[it]))/(r[it]*AA[it])
    end

    taxrev[1, it] = tauc[it]*CC[it]
    taxrev[2, it] = tauw[it]*w[it]*LL[it]
    taxrev[3, it] = taur[it]*r[it]*AA[it]
    taxrev[4, it] = sum(taxrev[1:3, it])

    # get budget balancing social security contribution
    taup[it] = PP[it]/(w[it]*LL[it])

end 


# computes the transition path of the economy
function get_transition()

    # initialize remaining variables
    if(!lsra_on)
        initialize_trn()
    else
        println("ITER    COMP_OLD  EFFICIENCY        DIFF")
    end

    # start timer
    #call tic()

    # iterate until value function converges
    for iter in 1:itermax

        # derive prices
        for it in 1:TT
            prices(it)
        end

        # solve the household problem
        for ij in JJ:-1:2
            solve_household(ij, 1)
        end
        for it in 1:TT
            solve_household(1, it)
        end

        # calculate the distribution of households over state space
        for it in 1:TT
            get_distribution(it)
        end

        # calculate lsra transfers if needed
        if(lsra_on)
            LSRA
        end

        # aggregate individual decisions
        for it in 1:TT
            aggregation(it)
        end

        # determine the government parameters
        for it in 1:TT
            government(it)
        end

        # get differences on goods markets
        itcheck = 0
        for it in 1:TT
            if(abs(DIFF[it]/YY[it])*100.0 < sig)
                itcheck = itcheck + 1
            end
        end

        #itmax = maxloc(abs(DIFF[1:TT]/YY[1:TT]), 1)
        itmax = argmax(abs.(DIFF[1:TT]./YY[1:TT]))

        # check for convergence and write screen output
        if(!lsra_on)

            check = iter > 1 && itcheck == TT && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0
            println(iter,"     ",round(digits = 5, HH[TT]),"   ", round(digits = 5, 5.0*KK[TT]/YY[TT]*100.0), "   ", round(digits = 5, CC[TT]/YY[TT]*100.0), "   ", round(digits = 5, II[TT]/YY[TT]*100.0), "   ", round(digits = 5, ((1.0+r[TT])^0.2-1.0)*100.0), "   ", round(digits = 5, w[TT]), "   ", round(digits = 8, DIFF[itmax]/YY[itmax]*100.0))

        else
            check = iter > 1 && itcheck == TT && lsra_comp/lsra_all > 0.99999 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0
            println(iter,"     ",round(digits = 5, lsra_comp/lsra_all*100.0), "    ", round(digits = 5, (Lstar^(1.0/(1.0-1.0/gamma))-1.0)*100.0), "     ", round(digits = 5, DIFF[itmax]/YY[itmax]*100.0))
        end

        # check for convergence
        if(check)
            for it in 1:TT
                if(!lsra_on)
                    output(it)
                end
            end
            #call output_summary()
            break
        end

    end

    #call toc
    for it in 1:TT
        if(!lsra_on)
            output(it)
        end
    end
    #call output_summary()

    #write(*,'(a/)')'ATTENTION: NO CONVERGENCE ###'

end 



# initializes transitional variables
function initialize_trn()

    println("TRANSITION PATH")

    println("ITER       H     K/Y     C/Y     I/Y       r       w        DIFF")

    #=
    # set up population structure
    for it in 1:TT
        pop[1, it] = (1.0+n_p)*pop[1, it-1]
        for ij in 2:JJ
            pop[ij, it] = pop[ij-1, it-1]
        end
    end

    for it in 1:TT
        for ij in 1:JJ
            m[ij, it] = pop[ij, it]/pop[1, it]
        end
    end
    =#
    
    for it in 1:TT

        taup[it] = taup[0]
        if (tax[it] == 1)
            tauc[it] = tauc[0]
        elseif (tax[it] == 2)
            tauw[it] = tauw[0]
            taur[it] = taur[0]
        elseif (tax[it] == 3)
            tauw[it] = tauw[0]
        else
            taur[it] = taur[0]
        end

        r[it] = r[0]
        rn[it] = r[it]*(1.0-taur[it])
        w[it] = w[0]
        wn[it] = w[it]*(1.0-tauw[it]-taup[it])
        p[it] = 1.0 + tauc[it]
        KK[it] = KK[0]
        AA[it] = AA[0]
        BB[it] = BB[0]
        LL[it] = LL[0]
        HH[it] = HH[0]
        YY[it] = YY[0]
        CC[it] = CC[0]
        II[it] = II[0]
        GG[it] = GG[0]
        INC[it] = INC[0]
        pen[:,it] = pen[:, 0]
        PP[it] = PP[0]
        taxrev[:,it] = taxrev[:, 0]
        c_coh[:, :, it] = c_coh[:, :, 0]
        l_coh[:, :, it] = l_coh[:, :, 0]
        y_coh[:, :, it] = y_coh[:, :, 0]
        a_coh[:, :, it] = a_coh[:, :, 0]
        aplus[:, :, :, :, :, it] = aplus[:, :, :, :, :, 0]
        c[:, :, :, :, :, it] = c[:, :, :, :, :, 0]
        l[:, :, :, :, :, it] = l[:, :, :, :, :, 0]
        phi[:, :, :, :, :, it] = phi[:, :, :, :, :, 0]
        VV[:, :, :, :, :, it] = VV[:, :, :, :, :, 0]
        RHS[:, :, :, :, :, it] = RHS[:, :, :, :, :, 0]
    end

end 
