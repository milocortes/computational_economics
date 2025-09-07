using DynamicProgrammingUtils
using Roots
using Printf
using Statistics

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

    # calculate tomorrows assets
    a_plus  = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)

    # get lsra transfer payment
    v_ind = v[ij_com, ia_com, ip_com, is_com, it_com]

    # calculate the wage rate
    wage = wn[it_com]*eff[ij_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+rn[it_com])*a[ia_com] + beq[ij_com, it_com] + pen[ij_com, it_com] + v_ind

    # determine labor
    if (ij_com < JR)
        lab_com = min( max( nu*(1.0-phi_l*Float64(ip_com)) + (1.0-nu)*(a_plus-available)/wage, 0.0) , 1.0-phi_l*Float64(ip_com)-1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max( (available + wage*lab_com - a_plus)/p[it_com] , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = max(varphi*RHS[ij_com+1, ial, ip_com, is_com, itp] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com, itp], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, ip_com, it_com)^(-gamma) - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, lab, im, it)

    return nu*(cons^nu*(1.0-lab-phi_l*Float64(im))^(1.0-nu))^(1.0-1.0/gamma)/(p[it]*cons)

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0),1.0-phi_l*Float64(ip)-1e-10)

    # get tomorrows year
    itp = year(it, ij, ij+1)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0

    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, is, itp] + (1.0-varphi)*EV[ij+1, iar, ip, is, itp], 1e-10)^(1.0-1.0/gamma)/(1.0-1.0/gamma)
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help-phi_l*Float64(ip))^(1.0-nu))^(1.0-1.0/gamma)/(1.0-1.0/gamma) + beta*psi[ij+1, itp]*valuefunc

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
    output(0)

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
    global pi_m 


    println("INITIAL EQUILIBRIUM")
    println("ITER     H     K/Y     C/Y     I/Y       r       w        DIFF")

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

    for it in 0:TT
        m[1, 0, it] = 0.44
        m[1, 1, it] = 1.0 - m[1, 0, it]
        
        itm = year2(it, -1)

        for ik in 0:NP
            for ij in 2:JJ
                m[ij, ik, it] = m[ij-1, ik, itm]*psi[ij, it]/(1.0+n_p)
            end
        end
    end

    for it in 0:TT

        GAM[1,it] = omega[1]
        itm = year2(it, -1)

        for ij in 2:JJ
            GAM[1,it] = GAM[1, it] + omega[ij]*sum(m[ij,:,it])
        end
        for ij in JJ:-1:1
            GAM[ij, it] = omega[ij]/GAM[1, it]
        end
    end

    # initialize asset grid
    grid_Cons_Grow(a,  NA+1, a_l, a_u, a_grow)

    # get initial guess for savings decision
    for ij in 1:JJ
        for ip in 0:NP
            for is in 1:NS
                @. aplus[ij, :, ip, is, 0] = max(a[:]/2.0, a[1]/2.0)
            end
        end
    end

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

    # initialize fixed effect
    dist_theta[0] = 0.44
    dist_theta[1] = 1 - dist_theta[0]

    theta[0]   = 1.0
    theta[1]   = exp(-varrho)
    #theta = exp.(theta)

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)

    # probability to have bad health when current health is good
    probability_bad_health_from_good_health_temp = zeros(JJ)
    grid_Cons_Equi(probability_bad_health_from_good_health_temp, JJ, 0.1, 0.4)
    pi_m[:, 0, 1] .= 0.05
    pi_m[:, 0, 0] .= 0.95

    # probability to have bad health when current health is bad
    probability_bad_health_from_bad_health_temp = zeros(JJ)

    grid_Cons_Equi(probability_bad_health_from_bad_health_temp, JJ, 0.6, 0.9)
    pi_m[:, 1, 1] .= 0.95
    pi_m[:, 1, 0] .= 0.05

    # tax and transfers
    tax   .= 2
    tauc  .= 0.275
    tauw  .= 0.0
    taur  .= 0.0
    taup  .= 0.1
    kappa .= 0.5
    gy    = 0.19
    by    = 0.60/5.0

    beq[:,0] .= 0.0
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

    # open files
    file_output = open("output.out", "w");
    file_summary = open("summary.out", "w");

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
    beq[JJ, it] = damp*GAM[JJ, it]*BQ[it] + (1.0-damp)*beq[JJ, it]#

    for ia in 0:NA
        for im in 0:NP
            aplus[JJ, ia, im, :, it] .= 0.0
            c[JJ, ia, im, :, it] .= ((1.0+rn[it])*a[ia] + beq[JJ, it] .+ pen[JJ, it] .+ v[JJ, ia, im, :, it])/p[it]
            l[JJ, ia, im, :, it] .= 0.0
            VV[JJ, ia, im, :, it] .= valuefunc(0.0, c[JJ, ia, im, 1, it],l[JJ, ia, im, 1, it], JJ, im, 1, it)
        end
    end

    # interpolate individual RHS
    interpolate(JJ, it)

    for ij in JJ-1:-1:ij_in

        it = year(it_in, ij_in, ij)

        #bequest for the olds
        beq[ij, it] = damp*GAM[ij, it]*BQ[it] + (1.0-damp)*beq[ij, it]#

        # check about how many is to iterate
        if (ij >= JR)
            ip_max = 1
            is_max = 1
        else
            ip_max = NP
            is_max = NS
        end

        for ia in 0:NA
            for ip in 0:NP
                # determine decision for zero assets at retirement without pension
                if (ij >= JR && ia == 0 && kappa[it] <= 1e-10)
                    aplus[ij, ia, ip, :, it] .= 0.0
                    c[ij, ia, ip, :, it] .= 0.0
                    l[ij, ia, ip, :, it] .= 0.0
                    VV[ij, ia, ip, :, it] .= valuefunc(0.0, 0.0, 0.0, ij, ip, 1, it)
                    continue
                end

                #for ip in 1:ip_max
                    for is in 1:is_max

                        # get initial guess for the individual choices
                        x_in = aplus[ij, ia, ip, is, it]

                        # set up communication variables
                        global ij_com = ij
                        global ia_com = ia
                        global ip_com = ip
                        global is_com = is
                        global it_com = it

                        # solve the household problem using rootfinding
                        #call fzero(x_in, foc, check)
                        x_root = fzero(foc, x_in)

                        # write screen output in case of a problem
                        #if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is, it

                        # check for borrowing constraint
                        if (x_root < 0.0)
                            x_root = 0.0
                            wage = wn[it]*eff[ij]*theta[ip]*eta[is]
                            v_ind = v[ij, ia, ip, is, it]
                            available = (1.0+rn[it])*a[ia] + beq[ij, it] + pen[ij, it] + v_ind
                            if (ij < JR)
                                global lab_com = min( max(nu*(1.0-phi_l*Float64(ip))-(1.0-nu)*available/wage , 0.0) , 1.0-phi_l*Float64(ip)-1e-10)
                            else
                                global lab_com = 0.0
                            end
                            global cons_com = max( (available + wage*lab_com)/p[it] , 1e-10)
                        end

                        # copy decisions
                        aplus[ij, ia, ip, is, it] = x_root
                        c[ij, ia, ip, is, it] = cons_com
                        l[ij, ia, ip, is, it] = lab_com
                        VV[ij, ia, ip, is, it] = valuefunc(x_root, cons_com, lab_com, ij, ip, is, it)

                    end

                    # copy decision in retirement age
                    if (ij >= JR)
                        aplus[ij, ia, ip, :, it] .= aplus[ij, ia, ip, 1, it]
                        c[ij, ia, ip, :, it] .= c[ij, ia, ip, 1, it]
                        l[ij, ia, ip, :, it] .= l[ij, ia, ip, 1, it]
                        VV[ij, ia, ip, :, it] .= VV[ij, ia, ip, 1, it]
                    end
                #end
            end
        end

        # interpolate individual RHS
        interpolate(ij, it)
        #println("Age :"*string(ij)*" DONE")
    end

end


# for calculating the rhs of the first order condition at age ij
function interpolate(ij, it)

    for ia in 0:NA
        for ip in 0:NP
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is, it] = 0.0
                EV[ij, ia, ip, is, it] = 0.0
                for ip_p in 0:NP
                    for is_p in 1:NS
                        chelp = max(c[ij, ia, ip_p, is_p, it],1e-10)
                        lhelp = max(l[ij, ia, ip_p, is_p, it],1e-10)
                        RHS[ij, ia, ip, is, it] = RHS[ij, ia, ip, is, it] + pi_m[ij, ip, ip_p]*pi[is, is_p]*margu(chelp, lhelp, ip_p, it)
                        EV[ij, ia, ip, is, it]  = EV[ij, ia, ip, is, it] + pi_m[ij, ip, ip_p]*pi[is, is_p]*VV[ij, ia, ip_p, is_p, it]
                    end
                end
                RHS[ij, ia, ip, is, it] = ((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ia, ip, is, it])^(-gamma)
                EV[ij, ia, ip, is, it] = ((1.0-1.0/gamma)*EV[ij, ia, ip, is, it])^(1.0/(1.0-1.0/gamma))
            end
        end
    end

end


# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, it] .= 0.0

    # get initial distribution in age 1
    for ip in 0:NP
        phi[1, 0, ip, is_initial, it] = dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 0:NP
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is, itm], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    ial = max(min(ial, NA-1), 0)
                    iar = max(min(iar, NA), 1)
                    varphi = max(min(varphi, 1.0), 0.0)

                    # redistribute households
                    for ip_p in 0:NP
                        for is_p in 1:NS
                            phi[ij, ial, ip_p, is_p, it] = phi[ij, ial, ip_p, is_p, it] + pi_m[ij, ip, ip_p]*pi[is, is_p]*varphi*phi[ij-1, ia, ip, is, itm]
                            phi[ij, iar, ip_p, is_p, it] = phi[ij, iar, ip_p, is_p, it] + pi_m[ij, ip, ip_p]*pi[is,is_p]*(1.0-varphi)*phi[ij-1,ia,ip,is,itm]
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

    m_coh = OffsetArray(zeros(JJ, NP+1), 1:JJ, 0:NP)

    # calculate fraction of good vs. bad health households
    for ij in 1:JJ
        for im in 0:NP
            frac_phi[ij, im, it] = sum(phi[ij, :, im, :, it])
        end
    end

    # calculate cohort aggregates
    c_coh[:, :, it]  .= 0.0
    l_coh[:, :, it]  .= 0.0
    y_coh[:, :, it]  .= 0.0
    a_coh[:, :, it]  .= 0.0
    VV_coh[:, :, it] .= 0.0
    m_coh[:, :]      .= 0.0
    FLC[:, :, it]     .= 0.0
    beq_coh[:, :, it] .= 0.0#

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 0:NP
                for is in 1:NS
                    c_coh[ij, ip, it] = c_coh[ij, ip, it] + c[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    l_coh[ij, ip, it] = l_coh[ij, ip, it] + l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    y_coh[ij, ip, it] = y_coh[ij, ip, it] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    a_coh[ij, ip, it] = a_coh[ij, ip, it] + a[ia]*phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]

                    # exclude households who die
                    if(ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        continue
                    end
                    if(aplus[ij, ia, ip, is, it] < 1e-4)
                        FLC[ij, ip, it] = FLC[ij, ip, it] + phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    end
                    VV_coh[ij, ip, it] = VV_coh[ij, ip, it] + VV[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    m_coh[ij, ip]      = m_coh[ij, ip] + phi[ij, ia, ip, is, it]/frac_phi[ij, ip,it]
                    beq_coh[ij, ip, it] = beq_coh[ij, ip, it] + a[ia]*(1.0 + rn[it])*(1.0 - psi[ij, it])*phi[ij, ia, ip, is, it]/psi[ij, it]/frac_phi[ij, ip,it] #
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

    frac_phi_adjusted = (m[:,:,it] ./ m[1,:,it]') .* frac_phi[:,:,it]

    m_adjusted[:,:, it] = frac_phi_adjusted

    for ij in 1:JJ
        for ip in 0:NP
            CC[it] = CC[it] + c_coh[ij, ip, it]*m_adjusted[ij, ip, it]
            LL[it] = LL[it] + y_coh[ij, ip, it]*m_adjusted[ij, ip, it]
            HH[it] = HH[it] + l_coh[ij, ip, it]*m_adjusted[ij, ip, it]
            AA[it] = AA[it] + a_coh[ij, ip, it]*m_adjusted[ij, ip, it]/psi[ij, it]#
            BQ[it] = BQ[it] + beq_coh[ij, ip, it]*m_adjusted[ij, ip, it]
            if (ij < JR)
                workpop = workpop + m_adjusted[ij, ip, it]
            end
        end
    end
    
    for it in 0:TT

        GAM[1,it] = omega[1]
        itm = year2(it, -1)

        for ij in 2:JJ
            GAM[1,it] = GAM[1, it] + omega[ij]*sum(m_adjusted[ij,:, it])
        end
        for ij in JJ:-1:1
            GAM[ij, it] = omega[ij]/GAM[1, it]
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
        for ip in 0:NP
            PP[it] = PP[it] + pen[ij, it]*m[ij, ip, it]
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
    global lsra_comp
    global lsra_all
    global Lstar
    global lsra_on

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
            LSRA()
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

            check = iter > 1 && itcheck == TT-1 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0

            fmt = "%4i " * "%8.2f "^6 * "%12.5f\n" #'(i4,6f8.2,f12.5)'

            @eval @printf $fmt $iter HH[TT] ([5.0*KK[TT], CC[TT], II[TT]]/YY[TT]*100)... ((1.0+r[TT])^0.2-1.0)*100.0 w[TT] DIFF[TT]/YY[TT]*100.0

        else

            check = iter > 1 && itcheck == TT-1 && lsra_comp/lsra_all > 0.99999 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0

            fmt = "%4i " * "%12.5f"^3 * "\n" # '(i4,3f12.5)'

            @eval @printf $fmt $iter  (lsra_comp/lsra_all*100.0) (Lstar^(1.0/(1.0-1.0/gamma))-1.0)*100.0 DIFF[$itmax]/YY[$itmax]*100.0
            
        end

        # check for convergence
        if(check)
            for it in 1:TT
                if(!lsra_on)
                    output(it)
                end
            end
            output_summary()
            break
        end

    end

    #call toc
    for it in 1:TT
        if(!lsra_on)
            output(it)
        end
    end
    output_summary()
    println("ATTENTION: NO CONVERGENCE ###")
    #write(*,'(a/)')'ATTENTION: NO CONVERGENCE ###'

end 



# initializes transitional variables
function initialize_trn()

    println("TRANSITION PATH")

    println("ITER       H     K/Y     C/Y     I/Y       r       w        DIFF")

    
    for it in 1:TT

        taup[it] = taup[0]
        tauc[it] = tauc[0]
        tauw[it] = tauw[0]
        taur[it] = taur[0]

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
        aplus[:, :, :, :, it] = aplus[:, :, :, :, 0]
        c[:, :, :, :, it] = c[:, :, :, :, 0]
        l[:, :, :, :, it] = l[:, :, :, :, 0]
        phi[:, :, :, :, it] = phi[:, :, :, :, 0]
        VV[:, :, :, :, it] = VV[:, :, :, :, 0]
        RHS[:, :, :, :, it] = RHS[:, :, :, :, 0]

        beq[:,it] = beq[:,0]
        BQ[it] = BQ[0]
    end

end 


# function for writing output
function output(it)

    global file_output
   # calculate cohort specific variances of logs
    exp_c = zeros(JJ)
    var_c = zeros(JJ)
    mas_c = zeros(JJ)
    exp_l = zeros(JJ)
    var_l = zeros(JJ)
    mas_l = zeros(JJ)
    exp_y = zeros(JJ)
    var_y = zeros(JJ)
    mas_y = zeros(JJ)


    for ij = 1:JJ
        for ia = 0:NA
            for ip = 0:NP
                for is = 1:NS

                    # consumption
                    if (c[ij, ia, ip, is, it] > 0.0)
                        temp = log(c[ij, ia, ip, is, it])
                        exp_c[ij] = exp_c[ij] + temp*phi[ij, ia, ip, is, it]
                        var_c[ij] = var_c[ij] + temp^2*phi[ij, ia, ip, is, it]
                        mas_c[ij] = mas_c[ij] + phi[ij, ia, ip, is, it]
                    end

                    if (l[ij, ia, ip, is, it] > 0.01)

                        # hours
                        temp = log(l[ij, ia, ip, is, it])
                        exp_l[ij] = exp_l[ij] + temp*phi[ij, ia, ip, is, it]
                        var_l[ij] = var_l[ij] + temp^2*phi[ij, ia, ip, is, it]
                        mas_l[ij] = mas_l[ij] + phi[ij, ia, ip, is, it]

                        # earnings
                        temp = log(w[it]*eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is, it])
                        exp_y[ij] = exp_y[ij] + temp*phi[ij, ia, ip, is, it]
                        var_y[ij] = var_y[ij] + temp^2*phi[ij, ia, ip, is, it]
                        mas_y[ij] = mas_y[ij] + phi[ij, ia, ip, is, it]
                    end
                end
            end
        end
    end


    exp_c = exp_c/max(maximum(mas_c), 1e-4)
    var_c = var_c/max(maximum(mas_c), 1e-4)
    exp_l = exp_l/max(maximum(mas_l), 1e-4)
    var_l = var_l/max(maximum(mas_l), 1e-4)
    exp_y = exp_y/max(maximum(mas_y), 1e-4)
    var_y = var_y/max(maximum(mas_y), 1e-4)
    var_c = var_c - exp_c.^2
    var_l = var_l - exp_l.^2
    var_y = var_y - exp_y.^2


    # Output
    @printf file_output "%s %3i \n\n" "EQUILIBRIUM YEAR " it 
    
    @printf file_output  "CAPITAL        K        A        B       BA        r     p.a.\n"
    
    fmt = " "^8*"%8.2f "^6*"\n"
    @eval @printf file_output  $fmt KK[$it] AA[$it] BB[$it] BA[$it]  r[$it]  ((1.0+r[$it])^(1.0/5.0)-1.0)*100.0

    fmt = "%s "*"%8.2f "^4 * "\n\n" 

    @eval @printf file_output $fmt  "(in %) " ([KK[$it], AA[$it], BB[$it], BA[$it]]/YY[$it]*500.0)...

    fmt = " "^8 * "%8.2f "^4 *"\n\n"
    @printf  file_output "%s \n" "LABOR           L       HH      INC        w"
    @eval @printf file_output  $fmt LL[$it] HH[$it]*100.0 INC[$it] w[$it]
    
    
    fmt_val = " "^8 * "%8.2f "^4 * "%8.3f\n"
    fmt_share = "%s " * "%8.2f "^4 * "%8.3f\n\n"
    @printf file_output "%s \n" "GOODS          Y       C       I       G    DIFF"
    @eval @printf file_output $fmt_val YY[$it] CC[$it] II[$it] GG[$it] DIFF[$it]
    @eval @printf file_output $fmt_share "(in %) " ([YY[$it], CC[$it], II[$it], GG[$it], DIFF[$it]]/YY[$it]*100.0)...

    fmt = " "^8 *"%8.2f"^6*"\n" 
    fmt_share = "%s"*"%8.2f"^6*"\n"
    fmt_rate = "%s"*"%8.2f"^3*"\n\n"
    @printf  file_output "%s \n" "GOV         TAUC    TAUW    TAUR    TOTAL      G       B"
    @eval @printf file_output $fmt taxrev[1:4, $it]... GG[$it] BB[$it]
    @eval @printf file_output $fmt_share  "(in %)  " (taxrev[1:4, $it]/YY[$it]*100)... ([GG[$it], BB[$it]*5.0]/YY[$it]*100.0)...
    @eval @printf file_output $fmt_rate  "(rate)  " ([tauc[$it], tauw[$it], taur[$it]]*100.0)...

    fmt_val = " "^8*"%8.2f"^3*"\n"
    fmt_rate = "%s "*"%8.2f"^3 * "\n\n"
    @printf  file_output "%s \n" "PENS        TAUP     PEN      PP"
    @eval @printf file_output $fmt_val taup[$it]*w[$it]*LL[$it] mean(pen[JR, $it]) PP[$it]
    @eval @printf file_output $fmt_rate "(in %) " ([taup[$it], kappa[$it], PP[$it]/YY[$it]]*100.0)... 

    fmt_val = " "^8*"%8.2f "^2*"\n"
    fmt_rate = "%s "*"%8.2f "^2*"\n\n"
    @printf file_output "%s \n" "LSRA          SV      BA"
    @eval @printf file_output $fmt_val SV[$it] BA[$it] 
    @eval @printf file_output $fmt_rate "(in %) " ([SV[$it], BA[$it]]/YY[$it]*100.0)...


    # check for the maximium grid point used
    iamax = check_grid(it)

    itp = year(it, 1, 2)

    @printf file_output "%s\n" "   IJ      CONS     LABOR  EARNINGS    INCOME    INCTAX      PENS    ASSETS    VAR(C)    VAR(L)    VAR(Y)      LSRA     VALUE     FLC       IAMAX     "

    fmt = "%3i "*"%10.3f "^12*"%10i %10i \n"

    for ij in 1:JJ
        for ip in 0:NP
            @eval @printf file_output $fmt $ij  sum(c_coh[$ij, $ip, $it])/INC[0]  sum(l_coh[$ij, $ip, $it]) ([w[$it]*sum(y_coh[$ij, $ip, $it]), wn[$it]*sum(y_coh[$ij, $ip, $it])+rn[$it]*sum(a_coh[$ij, $ip, $it]), tauw[$it]*w[$it]*sum(y_coh[$ij, $ip, $it])+taur[$it]*r[$it]*sum(a_coh[$ij, $ip, $it]), sum(pen[$ij, $it]-taup[$it]*w[$it]*y_coh[$ij, $ip, $it]) , 5.0*sum(a_coh[$ij, $ip, $it])]/INC[0])... $var_c[$ij] $var_l[$ij] $var_y[$ij] sum(v_coh[$ij, $ip, $it]) sum(VV_coh[$ij,$ip,$it]) sum(FLC[$ij, $ip, $it]) $iamax[$ij] 
        end
    end

    @printf file_output "%s \n\n" "--------------------------------------------------------------------"

end 



# subroutine that checks for the maximum gridpoint used
function check_grid(it)

    iamax = zeros(JJ)

    for ij = 1:JJ
        # check for the maximum asset grid point used at a certain age
        for ia = 0:NA
            for ip = 0:NP
                for is = 1:NS
                    if (phi[ij, ia, ip, is, it] > 1e-6)
                        iamax[ij] = ia
                    end
                end
            end     
        end
    end
    
    return iamax
end 


# writes summary output
function output_summary()

    HEV = OffsetArray(zeros( length(-(JJ-2):TT) ) , -(JJ-2):TT )
    mas = OffsetArray(zeros( length( -(JJ-2):0)), -(JJ-2):0)

    # aggregate ex post welfare changes of current generations
    HEV .= 0.0
    mas .= 0.0

    for ij = JJ:-1:2
        for ia = 0:NA
            for ip = 0:NP
                for is = 1:NS
                    if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        continue
                    end
                    HEV_help = ((VV[ij, ia, ip, is, 1]/max(VV[ij, ia, ip, is, 0], -1e10))^(1.0/egam)-1.0)*100.0
                    HEV[-(ij-2)] = HEV[-(ij-2)] + HEV_help*phi[ij, ia, ip, is, 1]
                    mas[-(ij-2)] = mas[-(ij-2)] + phi[ij, ia, ip, is, 1]
                end
            end
        end
    end

    HEV[-(JJ-2):0] = HEV[-(JJ-2):0]./parent(mas)

    # calculate ex ante welfare of future generations
    for it = 1:TT
        for ip = 0:NP
            HEV[it] = ((VV_coh[1, ip, it]/VV_coh[1, ip, 0])^(1.0/egam)-1.0)*100.0
        end
    end

    # headline
    @printf file_summary "%s \n"  "           A        K        L        H        r        w        C        I        Y        B       BA     tauc     tauw     taur     taup      HEV       DIFF"
    # current generations
    fmt = "%3i "*" "^135 * "%8.2f \n"
    for ij = -(JJ-2):-1
        @eval @printf file_summary $fmt $ij $HEV[$ij]
    end

    # future generations
    fmt = "%3i "*"%8.2f "^16 * "%10.5f \n"
    for it = 0:TT
        @eval @printf file_summary $fmt $it ([AA[$it]/AA[0]-1.0, KK[$it]/KK[0]-1.0, LL[$it]/LL[0]-1.0, HH[$it]-HH[0], (1.0+r[$it])^0.2-(1.0+r[0])^0.20, w[$it]/w[0]-1.0, CC[$it]/CC[0]-1.0, II[$it]/II[0]-1.0, YY[$it]/YY[0]-1.0, BB[$it]/BB[0]-1.0, BA[$it]/YY[$it], tauc[$it]-tauc[0], tauw[$it]-tauw[0], taur[$it]-taur[0], taup[$it]-taup[0]]*100.0)... $HEV[$it] DIFF[$it]/YY[$it]*100.0
    end

    if (lsra_on)
        @eval @printf file_summary "%s %12.6f \n" "EFFICIENCY GAIN: " (Lstar^(1.0/(1.0-1.0/gamma))-1.0)*100.0
    end

end 


# subroutine for calculating lsra payments
function LSRA()


    global lsra_comp
    global lsra_all
    global Lstar
    global lsra_on


    # initialize variables
    SV[:] .= 0.0
    v_coh[:, :, :] .= 0.0

    # initialize counters
    lsra_comp     = 0.0
    lsra_all      = 0.0

    for ij = 2:JJ
        for ia = 0:NA
            for ip = 0:NP
                for is = 1:NS

                    # do not do anything for an agent at retirement without pension and savings
                    if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        v[ij, ia, ip, is, 1] = 0.0
                        continue
                    end

                    # get today's utility
                    VV_1 = VV[ij, ia, ip, is, 1]

                    # get target utility
                    VV_0 = VV[ij, ia, ip, is, 0]

                    # get derivative of the value function
                    dVV_da = margu(c[ij, ia, ip, is, 1],l[ij, ia, ip, is, 1], ip, 1)*(1.0+rn[1])

                    # calculate change in transfers
                    v_tilde = (VV_0-VV_1)/dVV_da

                    # restrict z_tilde to income maximum
                    v_tilde = max(v_tilde, -((1.0+rn[1])*a[ia] + beq[ij,1] + pen[ij, 1] + wn[1]*eff[ij]*theta[ip]*eta[is]*0.99 + v[ij, ia, ip, is, 1]))
                    #v_tilde = max(v_tilde, -((1.0+rn[1])*a[ia] + beq[ij,1] + pen[ij, 1] + wn[1]*eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is, 1] + v[ij, ia, ip, is, 1]))
                    # check whether individual is already compensated
                    lsra_all = lsra_all + phi[ij, ia, ip, is, 1]/frac_phi[ij, ip,1]*m_adjusted[ij, ip, 1]
                    
                    if (abs((VV_1-VV_0)/VV_0)*100.0 < sig) 
                        lsra_comp = lsra_comp + phi[ij, ia, ip, is, 1]/frac_phi[ij, ip,1]*m_adjusted[ij, ip, 1]
                    end 

                    # calculate total transfer
                    v[ij, ia, ip, is, 1] = v[ij, ia, ip, is, 1] + damp*v_tilde

                    # aggregate transfers by cohort
                    v_coh[ij, ip, 1] = v_coh[ij, ip, 1] + v[ij, ia, ip, is, 1]*phi[ij, ia, ip, is, 1]/frac_phi[ij, ip,1]

                end
            end
        end
    end

    # aggregate transfers in year 1
    for ij = 2:JJ
        for ip = 0:NP
            SV[1] = SV[1] + v_coh[ij, ip, 1]*m_adjusted[ij, ip, 1]
        end
    end

    # initialize present value variables
    PV_t = 0.0
    PV_0 = 0.0
    PV_trans = 0.0

    # calculate present value of utility changes (in monetary values)
    for it = TT:-1:1
        for ip = 0:NP
            # get today's ex ante utility
            EVV_t = damp*VV_coh[1, ip, it]

            # get damped target utility
            EVV_0 = damp*VV_coh[1, ip, 0]

            # get derivative of expected utility function
            dEVV_da = 0.0

            for is = 1:NS
                dEVV_da = dEVV_da + margu(c[1, 0, ip, is, it], l[1, 0, ip, is, it], ip, it)*phi[1, 0, ip, is, it]/frac_phi[1, ip, it]*(1.0+rn[it])
            end

            # calculate present values
            if (it == TT)
                PV_t     = EVV_t/dEVV_da    *(1.0+r[it])/(r[it]-n_p)
                PV_0     = EVV_0/dEVV_da    *(1.0+r[it])/(r[it]-n_p)
                PV_trans = v[1, 0, ip, 1, it]*(1.0+r[it])/(r[it]-n_p)
            else
                PV_t     = PV_t    *(1.0+n_p)/(1.0+r[it+1]) + EVV_t/dEVV_da
                PV_0     = PV_0    *(1.0+n_p)/(1.0+r[it+1]) + EVV_0/dEVV_da
                PV_trans = PV_trans*(1.0+n_p)/(1.0+r[it+1]) + v[1, 0, ip, 1, it]
            end
        end
    end

    # calculate the constant utility gain/loss for future generations
    Lstar = (PV_t-PV_trans-SV[1])/PV_0

    # calculate compensation payments for future cohorts
    for it = TT:-1:1
        for ip = 0:NP
            # get today's ex ante utility
            EVV_t = damp*VV_coh[1, ip, it]

            # get target utility
            EVV_0 = damp*VV_coh[1, ip, 0]*Lstar

            # get derivative of expected utility function
            dEVV_da = 0.0

            for is = 1:NS
                dEVV_da = dEVV_da + margu(c[1, 0, ip, is, it], l[1, 0, ip, is, it], ip, it)*phi[1, 0, ip, is, it]/frac_phi[1, ip, it]*(1.0+rn[it])
            end
        
            # compute change in transfers (restricted)
            v_tilde = (EVV_0-EVV_t)/dEVV_da

            # calculate cohort transfer level
            #v[1, 0, 0, :, it] = v[1, 0, 0, :, it] .+ v_tilde
            v[1, 0, ip, :, it] = v[1, 0, ip, :, it] .+ v_tilde

            # aggregate transfers
            #v_coh[1, ip , it] = v[1, 0, 1, 1, it]
            v_coh[1, ip , it] = v[1, 0, ip, 1, it]
            SV[it] = SV[it] + v_coh[1, ip, it]*m_adjusted[1, ip, it]
        end
    end

    # determine sequence of LSRA debt/savings
    BA[2] = SV[1]/(1.0+n_p)
    for it = 3:TT
        BA[it] = ((1.0+r[it-1])*BA[it-1] + SV[it-1])/(1.0+n_p)
    end

end 

