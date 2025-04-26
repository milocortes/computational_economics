using DynamicProgrammingUtils
using Roots



# the first order condition
function foc(x_in)
    global ij_com
    global is_com
    global cons_com
    global lab_com

    # calculate tomorrows assets
    a_plus  = x_in

    # calculate the wage rate
    wage = wn*eff[ij_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+rn)*a[ia_com] + pen[ij_com]

    # determine labor
    if (ij_com < JR)
        global lab_com = min( max( nu + (1.0-nu)*(a_plus - available)/wage, 0.0), 1.0-1e-10)
    else
        global lab_com = 0.0
    end

    # calculate consumption
    global cons_com = max( (available + wage*lab_com - a_plus)/p , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, ip_com, is_com] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com]

    # calculate first order condition for consumption
    return margu(cons_com, lab_com)^(-gamma) - tomorrow

end 


# calculates marginal utility of consumption
function margu(cons, lab)

    return nu/p*(cons^nu*(1.0-lab)^(1.0-nu))^egam/cons

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is)


    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0),1.0-1e-10)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, is] + (1.0-varphi)*EV[ij+1, iar, ip, is], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^egam/egam + beta*valuefunc

end 


# computes the initial steady state of the economy
function get_SteadyState()


    # iterate until value function converges
    for iter in 1:itermax

        # derive prices
        prices()

        # solve the household problem
        solve_household()

        # calculate the distribution of households over state space
        get_distribution()

        # aggregate individual decisions over cohorts
        aggregation()

        # determine the government parameters
        government()

        println(iter,"     ",round(digits = 2, 5.0*KK/YY*100),"   ", round(digits = 2, CC/YY*100.0),"   ", round(digits = 2, II/YY*100.0), "   ", round(digits = 2, r),"   ", round(digits = 2, w), "   ", round(digits = 2, DIFF/YY*100.0))


        if (abs(DIFF/YY)*100.0 < sig)
            #call toc
            #call output()
            return
        end
    end
end 


# initializes the remaining model parameters and variables
function initialize()
    global a 
    global aplus
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
    global KK 
    global LL 
    global YY
    global II
    global GG
    global BB
    global pen
    
    println("INITIAL EQUILIBRIUM")
    println("ITER     K/Y     C/Y     I/Y       r       w        DIFF")

    # set up population structure
    for ij in 1:JJ
        m[ij] = (1.0+n_p)^(1.0-ij)
    end

    # initialize asset grid
    grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)

    # get initial guess for savings decision
    for ij in 1:JJ
        for ip in 1:NP
            for is in 1:NS
                aplus[ij, :, ip, is] .= maximum(a[:]/2.0)
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
    dist_theta .= 1.0/Float64(NP)
    theta[1]   = -sqrt(sigma_theta)
    theta[2]   = sqrt(sigma_theta)
    theta = exp.(theta)

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
    eta = exp.(eta)

    # tax and transfers
    tax   = 2
    tauc  = 0.075
    tauw  = 0.0
    taur  = 0.0
    taup  = 0.1
    kappa = 0.5
    gy    = 0.19
    by    = 0.60/5.0

    # initial guesses for macro variables
    KK = 1.0
    LL = 1.0
    YY = 1.0
    II = (n_p+delta)*KK

    GG = gy*YY
    BB = by*YY

    pen .= 0.0
    pen[JR:JJ] .= kappa

end 


# subroutine for calculating prices
function prices()
    global r 
    global w 
    global rn 
    global wn 
    global p 

    r = Omega*alpha*(KK/LL)^(alpha-1.0)-delta
    w = Omega*(1.0-alpha)*(KK/LL)^alpha
    rn = r*(1.0-taur)
    wn = w*(1.0-tauw-taup)
    p = 1.0 + tauc

end 


# determines the solution to the household optimization problem
function solve_household()
    global cons_com
    global lab_com


    # get decision in the last period of life
    for ia in 0:NA
        aplus[JJ, ia, :, :] .= 0.0
        c[JJ, ia, :, :] .= ((1.0+rn)*a[ia] + pen[JJ])/p
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
            if (ij >= JR && ia == 0 && kappa <= 1e-10)
                aplus[ij, ia, :, :] .= 0.0
                c[ij, ia, :, :] .= 0.0
                l[ij, ia, :, :] .= 0.0
                V[ij, ia, :, :] .= valuefunc(0.0, 0.0, 0.0, ij, 1, 1)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij, ia, ip, is]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is

                    # solve the household problem using rootfinding
                    #call fzero(x_in, foc, check)
                    x_root = fzero(foc, x_in)

                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is

                    # check for borrowing constraint
                    if (x_root < 0.0)
                        x_root = 0.0
                        wage = wn*eff[ij]*theta[ip]*eta[is]
                        available = (1.0+rn)*a[ia] + pen[ij]
                        if (ij < JR)
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        global cons_com = max( (available + wage*lab_com)/p , 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, ip, is] = x_root
                    c[ij, ia, ip, is] = cons_com
                    l[ij, ia, ip, is] = lab_com
                    V[ij, ia, ip, is] = valuefunc(x_root, cons_com, lab_com, ij, ip, is)

                end

                # copy decision in retirement age
                if (ij >= JR)
                    aplus[ij, ia, :, :] .= aplus[ij, ia, 1, 1]
                    c[ij, ia, :, :] .= c[ij, ia, 1, 1]
                    l[ij, ia, :, :] .= l[ij, ia, 1, 1]
                    V[ij, ia, :, :] .= V[ij, ia, 1, 1]
                end
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
                    chelp = max(c[ij, ia, ip, is_p],1e-10)
                    lhelp = max(l[ij, ia, ip, is_p],1e-10)
                    RHS[ij, ia, ip, is] = RHS[ij, ia, ip, is] + pi[is, is_p]*margu(chelp, lhelp)
                    EV[ij, ia, ip, is]  = EV[ij, ia, ip, is] + pi[is, is_p]*V[ij, ia, ip, is_p]
                end
                RHS[ij, ia, ip, is] = ((1.0+rn)*beta*RHS[ij, ia, ip, is])^(-gamma)
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
    for ij = 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    ial = max(min(ial, NA-1), 0)
                    iar = max(min(iar, NA), 1)
                    varphi = max(min(varphi, 1,0), 0.0)

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
    global LL_old = LL

    # calculate cohort aggregates
    c_coh[:] .= 0.0
    l_coh[:] .= 0.0
    y_coh[:] .= 0.0
    a_coh[:] .= 0.0
    v_coh[:] .= 0.0

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    c_coh[ij] = c_coh[ij] + c[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    l_coh[ij] = l_coh[ij] + l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    y_coh[ij] = y_coh[ij] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    a_coh[ij] = a_coh[ij] + a[ia]*phi[ij, ia, ip, is]
                    v_coh[ij] = v_coh[ij] + V[ij, ia, ip, is]*phi[ij, ia, ip, is]
                end
            end
        end
    end

    # calculate aggregate quantities
    global CC = 0.0
    global LL = 0.0
    global HH = 0.0
    global AA = 0.0
    global workpop = 0.0

    for ij in 1:JJ
        global CC = CC + c_coh[ij]*m[ij]
        global LL = LL + y_coh[ij]*m[ij]
        global HH = HH + l_coh[ij]*m[ij]
        global AA = AA + a_coh[ij]*m[ij]
        if (ij < JR)
            global workpop = workpop + m[ij]
        end
    end

    # damping and other quantities
    global KK = damp*(AA-BB) + (1.0-damp)*KK
    global LL = damp*LL + (1.0-damp)*LL_old
    global II = (n_p+delta)*KK
    global YY = Omega * KK^alpha * LL^(1.0-alpha)

    # get average income and average working hours
    global INC = w*LL/workpop
    global HH  = HH/workpop

    # get difference on goods market
    global DIFF = YY-CC-II-GG

    
end 



# subroutine for calculating government parameters
function government()

    # set government quantities and pension payments
    if(!reform_on)
        global GG = gy*YY
        global BB = by*YY
    end

    # calculate government expenditure
    global expend = GG + (1.0+r)*BB - (1.0+n_p)*BB

    # get budget balancing tax rate
    if (tax == 1)
        global tauc = (expend - (tauw*w*LL + taur*r*AA))/CC
        global p    = 1.0 + tauc
    elseif (tax == 2)
        global tauw = (expend - tauc*CC)/(w*LL + r*AA)
        global taur = tauw
    elseif (tax == 3)
        global tauw = (expend - (tauc*CC + taur*r*AA))/(w*LL)
    else
        global taur = (expend - (tauc*CC + tauw*w*LL))/(r*AA)
    end

    taxrev[1] = tauc*CC
    taxrev[2] = tauw*w*LL
    taxrev[3] = taur*r*AA
    taxrev[4] = sum(taxrev[1:3])

    # get budget balancing social security contribution
    pen[JR:JJ] .= kappa*INC
    global PP = 0.0
    for ij in JR:JJ
        global PP = PP + pen[ij]*m[ij]
    end

    global taup = PP/(w*LL)

end 