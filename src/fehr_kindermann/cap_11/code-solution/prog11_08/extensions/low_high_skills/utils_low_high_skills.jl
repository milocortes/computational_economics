using DynamicProgrammingUtils
using Roots
using Printf
using Statistics

# calculates year at which age ij agent is ij_p
function year(it, ij, ijp)

    year = it + ijp - ij

    if (it == 0 || year <= 0)
        year = 0
    end

    if (it == TT || year >= TT)
        year = TT
    end

    return year
end

# function which computes the year in which the household lives
function year2(it, addit)


    year2 = it + addit

    if (year2 > TT)
        year2 = TT
    end

    if (year2 < 0)
        year2 = 0
    end

    if (it == 0)
        year2 = 0
    end
    if (it == TT)
        year2 = TT
    end

    return year2
end

function foc(x_in)
    global ij_com
    global ia_com
    global ip_com
    global is_com
    global it_com
    global cons_com
    global lab_com
    global ir_com
    global epplus_com
    global ik_com

    # calculate tomorrows assets
    a_plus = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com + 1)

    # get lsra transfer payment
    v_ind = v[ij_com, ik_com, ia_com, ir_com, ip_com, is_com, it_com]

    # calculate the marginal wage rate
    wage = w[it_com] * eff[ij_com, ik_com] * theta[ip_com] * eta[is_com]
    wagen = wn[it_com] * eff[ij_com, ik_com] * theta[ip_com] * eta[is_com]

    # calculate available resources
    available = (1.0 + rn[it_com]) * a[ia_com] + penp[ij_com, ik_com, it_com, ir_com] + v_ind

    # determine labor
    if (ij_com < JR)
        lab_com = min(max(nu + (1.0 - nu) * (a_plus - available) / (wage * (1.0 - tauw[it_com] - tau_impl[ij_com, it_com])), 0.0), 1.0 - 1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max((available + wagen * lab_com - a_plus) / p[it_com], 1e-10)

    # pension system
    if (ij_com >= JR)
        epplus_com = ep[ir_com]
    else
        epplus_com = (ep[ir_com] * Float64(ij_com - 1) + lambda[it_com] + (1.0 - lambda[it_com]) * wage * lab_com / INC[it_com]) / Float64(ij_com)
    end


    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi_a = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    irl, irr, varphi_r = linint_Grow(epplus_com, ep_l, ep_u, ep_grow, NR, ial_ep, iar_ep, varphi_ep)


    tomorrow = max(varphi_a * varphi_r * RHS[ij_com+1, ik_com, ial, irl, ip_com, is_com, itp] + varphi_a * (1.0 - varphi_r) * RHS[ij_com+1, ik_com, ial, irr, ip_com, is_com, itp] + (1.0 - varphi_a) * varphi_r * RHS[ij_com+1, ik_com, iar, irl, ip_com, is_com, itp] + (1.0 - varphi_a) * (1.0 - varphi_r) * RHS[ij_com+1, ik_com, iar, irr, ip_com, is_com, itp], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, it_com)^(-gamma) - tomorrow

end


# calculates marginal utility of consumption
function margu(cons, lab, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0), 1.0 - 1e-10)

    return nu * (c_help^nu * (1.0 - l_help)^(1.0 - nu))^egam / (p[it] * c_help)

end

# calculates the value function
function valuefunc(a_plus, ep_plus, cons, lab, ij, ik, ip, is, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0), 1.0 - 1e-10)

    # get tomorrows year
    itp = year(it, ij, ij + 1)

    # get tomorrows utility
    ial, iar, varphi_a = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)
    irl, irr, varphi_r = linint_Grow(ep_plus, ep_l, ep_u, ep_grow, NR, ial_ep, iar_ep, varphi_ep)


    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if (ij < JJ)
        valuefunc = max(varphi_a * varphi_r * EV[ij+1, ik, ial, irl, ip, is, itp] + varphi_a * (1.0 - varphi_r) * EV[ij+1, ik, ial, irr, ip, is, itp] + (1.0 - varphi_a) * varphi_r * EV[ij+1, ik, iar, irl, ip, is, itp] + (1.0 - varphi_a) * (1.0 - varphi_r) * EV[ij+1, ik, iar, irr, ip, is, itp], 1e-10)^egam / egam
    end

    # add todays part and discount
    return (c_help^nu * (1.0 - l_help)^(1.0 - nu))^egam / egam + beta * valuefunc

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

        # get implicit tax rates
        implicit_taxes(0)

        # solve the household problem
        solve_household(1, 0)

        # calculate the distribution of households over state space
        get_distribution(0)

        # aggregate individual decisions
        aggregation(0)

        # determine the government parameters
        government(0)

        fmt = "%4i " * "%8.2f "^6 * "%12.5f\n"#'(i4,6f8.2,f12.5)'
        @eval @printf $fmt $iter HH[0] ([5.0*KK[0], CC[0], II[0]]/YY[0]*100)... ((1.0+r[0])^0.2-1.0)*100.0 w[0] DIFF[0]/YY[0]*100.0
        #println(iter,"     ",round(digits = 2, HH[0]),"   ", round(digits = 2, 5.0*KK[0]/YY[0]*100.0), "   ", round(digits = 2, CC[0]/YY[0]*100.0), "   ", round(digits = 2, II[0]/YY[0]*100.0), "   ", round(digits = 2, ((1.0+r[0])^0.2-1.0)*100.0), "   ", round(digits = 2, w[0]), "   ", round(digits = 6, DIFF[0]/YY[0]*100.0))

        if (abs(DIFF[0] / YY[0]) * 100.0 < sig)
            break
        end
    end

    #call toc
    #call output(0)

    #write(*,'(a/)')'ATTENTION: NO CONVERGENCE ###'

end



# initializes the remaining model parameters and variables
function initialize()

    println("INITIAL EQUILIBRIUM")
    println("ITER     H     K/Y     C/Y     I/Y       r       w        DIFF")
    ## Population parameters
    global m, pop

    ## Asset arrays
    global a, a_plus

    ## Shock process parameters
    global eff, dist_theta, theta, pi, eta

    ## Taxes
    global tax, tauc, tauw, taur, taup

    ## Pensions
    global kappa, pen, ep, penp, lambda

    ## Goverment parameters
    global gy, by

    ## Aggregated parameters
    global KK, LL, YY, II, GG, BB

    # set up population structure
    #=
    for ij in 1:JJ
        pop[ij, 0] = 1.0 / (1.0 + n_p)^(ij - 1)
    end

    for ij in 1:JJ
        m[ij, 0] = pop[ij, 0] / pop[1, 0]
    end
    =#
    for it in 0:TT
        m[1, 1, it] = 1 - 0.6
        m[1, 2, it] = 1.0 - m[1, 1, it]

        itm = year2(it, -1)

        for ik in 1:SS
            for ij in 2:JJ
                m[ij, ik, it] = m[ij-1, ik, itm] / (1.0 + n_p)
            end
        end
    end

    # initialize asset and pension grid
    grid_Cons_Grow(a, NA + 1, a_l, a_u, a_grow)
    grid_Cons_Grow(ep, NR + 1, ep_l, ep_u, ep_grow)

    # get initial guess for savings decision and labor supply
    for ij in 1:JJ
        for ik in 1:SS
            for ir in 0:NR
                for ip in 1:NP
                    for is in 1:NS
                        @. aplus[ij, ik, :, ir, ip, is, 0] = max(a[:] / 2.0, a[1] / 2.0) # max(a[:]/2.0, a[1]/2.0)
                        @. c[ij, ik, :, ir, ip, is, 0] = max(a[:] / 2.0, a[1] / 2.0)
                    end
                end
            end
        end
    end

    # initialize age earnings process
    eff[1, 2] = 1.0000
    eff[2, 2] = 1.3527
    eff[3, 2] = 1.6952
    eff[4, 2] = 1.8279
    eff[5, 2] = 1.9606
    eff[6, 2] = 1.9692
    eff[7, 2] = 1.9692
    eff[8, 2] = 1.9392
    eff[9, 2] = 1.9007
    eff[JR:JJ, 2] .= 0.0

    eff[:, 1] = eff[:, 2] * 0.8

    # initialize fixed effect
    dist_theta .= 1.0 / Float64(NP)
    theta[1] = -sqrt(sigma_theta)
    theta[2] = sqrt(sigma_theta)
    theta .= exp.(theta)

    # calculate the shock process
    pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0)
    eta = exp.(eta)

    # tax and transfers
    tax .= 2
    tauc .= 0.075
    tauw .= 0.2088
    taur .= 0.2088
    taup .= 0.12
    kappa .= 0.50
    #lambda = 1d0
    lambda .= 0.0
    penp .= 0.0
    gy = 0.19
    by = 0.60 / 5.0

    # initial guesses for macro variables
    KK .= 5.14
    LL .= 5.34
    YY .= 8.42
    II .= (n_p + delta) * KK
    INC .= 0.72

    GG .= gy * YY[0]
    BB .= by * YY[0]
    PP .= 0.0

    # open files
    file_output = open("output.out", "w");
    file_summary = open("summary.out", "w");
end


# function for calculating prices
function prices(it)

    r[it] = Omega * alpha * (KK[it] / LL[it])^(alpha - 1.0) - delta
    w[it] = Omega * (1.0 - alpha) * (KK[it] / LL[it])^alpha
    rn[it] = r[it] * (1.0 - taur[it])
    wn[it] = w[it] * (1.0 - tauw[it] - taup[it])
    p[it] = 1.0 + tauc[it]

end


# determines the solution to the household optimization problem
function solve_household(ij_in, it_in)
    global cons_com
    global lab_com

    # get decision in the last period of life
    it = year(it_in, ij_in, JJ)

    for ik in 1:SS
        for ia in 0:NA
            aplus[JJ, ik, ia, :, :, :, it] .= 0.0
            epplus[JJ, ik, ia, :, :, :, it] .= 0.0
            l[JJ, ik, ia, :, :, :, it] .= 0.0

            for ir in 0:NR
                penp[JJ, ik, it, ir] = kappa[it] * INC[it] * ep[ir]
                @. c[JJ, ik, ia, ir, :, :, it] .= ((1.0 + rn[it]) * a[ia] + penp[JJ, ik, it, ir] + v[JJ, ik, ia, ir, :, :, it]) / p[it]
                VV[JJ, ik, ia, ir, :, :, it] .= valuefunc(0.0, 0.0, c[JJ, ik, ia, ir, 1, 1, it], l[JJ, ik, ia, ir, 1, 1, it], JJ, ik, 1, 1, it)
            end
        end
    end

    # interpolate individual expectations
    interpolate(JJ, 1, it)
    interpolate(JJ, 2, it)

    for ij in JJ-1:-1:ij_in
        #println(ij)
        for ik in 1:SS
            it = year(it_in, ij_in, ij)

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
                    aplus[ij, ik, ia, :, :, :, it] .= 0.0
                    epplus[ij, ik, ia, :, :, :, it] .= 0.0
                    c[ij, ik, ia, :, :, :, it] .= 0.0
                    l[ij, ik, ia, :, :, :, it] .= 0.0
                    VV[ij, ik, ia, :, :, :, it] .= valuefunc(0.0, 0.0, 0.0, 0.0, ij, ik, 1, 1, it)
                    continue
                end

                for ir in 0:NR
                    for ip in 1:ip_max
                        for is in 1:is_max
                            # pension system
                            if (ij >= JR)
                                penp[ij, ik, it, ir] = kappa[it] * INC[it] * ep[ir]
                            else
                                penp[ij, ik, it, ir] = 0.0
                            end

                            # get initial guess for the individual choices
                            x_in = aplus[ij, ik, ia, ir, ip, is, it]

                            # set up communication variables
                            global ij_com = ij
                            global ia_com = ia
                            global ir_com = ir
                            global ip_com = ip
                            global is_com = is
                            global it_com = it
                            global ik_com = ik

                            # solve the household problem using rootfinding
                            #fzero(x_in, foc, check)
                            x_root = fzero(foc, x_in)

                            # write screen output in case of a problem
                            #if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ia, ir, ip, is, it

                            # check for borrowing constraint
                            if (x_root < 0.0)
                                x_root = 0.0
                                temp = foc(x_root)
                            end

                            # copy decisions
                            aplus[ij, ik, ia, ir, ip, is, it] = x_root
                            epplus[ij, ik, ia, ir, ip, is, it] = epplus_com
                            c[ij, ik, ia, ir, ip, is, it] = cons_com
                            l[ij, ik, ia, ir, ip, is, it] = lab_com
                            VV[ij, ik, ia, ir, ip, is, it] = valuefunc(x_root, epplus_com, cons_com, lab_com, ij, ik, ip, is, it)
                        end

                        # copy decision in retirement age
                        if (ij >= JR)
                            aplus[ij, ik, ia, ir, :, :, it] .= aplus[ij, ik, ia, ir, 1, 1, it]
                            epplus[ij, ik, ia, ir, :, :, it] .= epplus[ij, ik, ia, ir, 1, 1, it]
                            c[ij, ik, ia, ir, :, :, it] .= c[ij, ik, ia, ir, 1, 1, it]
                            l[ij, ik, ia, ir, :, :, it] .= l[ij, ik, ia, ir, 1, 1, it]
                            VV[ij, ik, ia, ir, :, :, it] .= VV[ij, ik, ia, ir, 1, 1, it]
                        end
                    end
                end
            end

            # interpolate individual expectations
            interpolate(ij, ik, it)
        end
    end

end


# for interpolating individual expectations
function interpolate(ij, ik, it)

    for ia in 0:NA
        for ir in 0:NR
            for ip in 1:NP
                for is in 1:NS

                    # calculate RHS and EV
                    RHS[ij, ik, ia, ir, ip, is, it] = 0.0
                    EV[ij, ik, ia, ir, ip, is, it] = 0.0

                    for is_p in 1:NS
                        chelp = max(c[ij, ik, ia, ir, ip, is_p, it], 1e-10)
                        lhelp = max(l[ij, ik, ia, ir, ip, is_p, it], 1e-10)
                        RHS[ij, ik, ia, ir, ip, is, it] = RHS[ij, ik, ia, ir, ip, is, it] + pi[is, is_p] * margu(chelp, lhelp, it)
                        EV[ij, ik, ia, ir, ip, is, it] = EV[ij, ik, ia, ir, ip, is, it] + pi[is, is_p] * VV[ij, ik, ia, ir, ip, is_p, it]
                    end
                    RHS[ij, ik, ia, ir, ip, is, it] = ((1.0 + rn[it]) * beta * RHS[ij, ik, ia, ir, ip, is, it])^(-gamma)
                    EV[ij, ik, ia, ir, ip, is, it] = (egam * EV[ij, ik, ia, ir, ip, is, it])^(1.0 / egam)
                end
            end
        end
    end

end


# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, :, :, it] .= 0.0

    # get initial distribution in age 1
    for ik in 1:SS
        for ip in 1:NP
            phi[1, ik, 0, 0, ip, is_initial, it] = dist_theta[ip]
        end
    end

    # successively compute distribution over ages
    for ij in 2:JJ
        for ik in 1:SS
            # iterate over yesterdays gridpoints
            for ia in 0:NA
                for ir in 0:NR
                    for ip in 1:NP
                        for is in 1:NS

                            # interpolate yesterday's savings decision
                            ial, iar, varphi_a = linint_Grow(aplus[ij-1, ik, ia, ir, ip, is, itm], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)
                            irl, irr, varphi_r = linint_Grow(epplus[ij-1, ik, ia, ir, ip, is, itm], ep_l, ep_u, ep_grow, NR, ial_ep, iar_ep, varphi_ep)

                            # restrict values to grid just in case
                            ial = min(ial, NA)
                            iar = min(iar, NA)
                            varphi_a = max(0.0, min(varphi_a, 1.0))
                            irl = min(irl, NR)
                            irr = min(irr, NR)
                            varphi_r = max(0.0, min(varphi_r, 1.0))

                            # redistribute households
                            for is_p in 1:NS
                                phi[ij, ik, ial, irl, ip, is_p, it] = phi[ij, ik, ial, irl, ip, is_p, it] + pi[is, is_p] * varphi_a * varphi_r * phi[ij-1, ik, ia, ir, ip, is, itm]
                                phi[ij, ik, ial, irr, ip, is_p, it] = phi[ij, ik, ial, irr, ip, is_p, it] + pi[is, is_p] * varphi_a * (1.0 - varphi_r) * phi[ij-1, ik, ia, ir, ip, is, itm]
                                phi[ij, ik, iar, irl, ip, is_p, it] = phi[ij, ik, iar, irl, ip, is_p, it] + pi[is, is_p] * (1.0 - varphi_a) * varphi_r * phi[ij-1, ik, ia, ir, ip, is, itm]
                                phi[ij, ik, iar, irr, ip, is_p, it] = phi[ij, ik, iar, irr, ip, is_p, it] + pi[is, is_p] * (1.0 - varphi_a) * (1.0 - varphi_r) * phi[ij-1, ik, ia, ir, ip, is, itm]
                            end
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
    c_coh[:, :, it] .= 0.0
    l_coh[:, :, it] .= 0.0
    y_coh[:, :, it] .= 0.0
    a_coh[:, :, it] .= 0.0
    pen[:, :, it] .= 0.0
    VV_coh[:, :, it] .= 0.0
    m_coh[:, :] .= 0.0
    FLC[:, :, it] .= 0.0

    for ij in 1:JJ
        for ik in 1:SS
            for ia in 0:NA
                for ir in 0:NR
                    for ip in 1:NP
                        for is in 1:NS
                            c_coh[ij, ik, it] = c_coh[ij, ik, it] + c[ij, ik, ia, ir, ip, is, it] * phi[ij, ik, ia, ir, ip, is, it]
                            l_coh[ij, ik, it] = l_coh[ij, ik, it] + l[ij, ik, ia, ir, ip, is, it] * phi[ij, ik, ia, ir, ip, is, it]
                            y_coh[ij, ik, it] = y_coh[ij, ik, it] + eff[ij, ik] * theta[ip] * eta[is] * l[ij, ik, ia, ir, ip, is, it] * phi[ij, ik, ia, ir, ip, is, it]
                            a_coh[ij, ik, it] = a_coh[ij, ik, it] + a[ia] * phi[ij, ik, ia, ir, ip, is, it]
                            pen[ij, ik, it] = pen[ij, ik, it] + penp[ij, ik, it, ir] * phi[ij, ik, ia, ir, ip, is, it]
                            # exclude households who dies
                            if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                                continue
                            end
                            if (aplus[ij, ik, ia, ir, ip, is, it] < 1e-4)
                                FLC[ij, ik, it] = FLC[ij, ik, it] + phi[ij, ik, ia, ir, ip, is, it]
                            end
                            VV_coh[ij, ik, it] = VV_coh[ij, ik, it] + VV[ij, ik, ia, ir, ip, is, it] * phi[ij, ik, ia, ir, ip, is, it]
                            m_coh[ij, ik] = m_coh[ij, ik] + phi[ij, ik, ia, ir, ip, is, it]
                        end
                    end
                end
            end
        end
    end

    # normalize VV_coh (because hh excluded)
    VV_coh[:, :, it] = VV_coh[:, :, it] ./ m_coh
    FLC[:, :, it] = FLC[:, :, it] ./ m_coh

    # calculate aggregate quantities
    CC[it] = 0.0
    LL[it] = 0.0
    HH[it] = 0.0
    AA[it] = 0.0
    workpop[it] = 0.0

    for ij in 1:JJ
        for ik in 1:SS
            CC[it] = CC[it] + c_coh[ij, ik, it] * m[ij, ik, it]
            LL[it] = LL[it] + y_coh[ij, ik, it] * m[ij, ik, it]
            HH[it] = HH[it] + l_coh[ij, ik, it] * m[ij, ik, it]
            AA[it] = AA[it] + a_coh[ij, ik, it] * m[ij, ik, it]
            if (ij < JR)
                workpop[it] = workpop[it] + m[ij, ik, it]
            end
        end
    end

    # damping and other quantities
    KK[it] = damp * (AA[it] - BB[it] - BA[it]) + (1.0 - damp) * KK[it]
    LL[it] = damp * LL[it] + (1.0 - damp) * LL_old
    II[it] = (1.0 + n_p) * KK[itp] - (1.0 - delta) * KK[it]
    YY[it] = Omega * KK[it]^alpha * LL[it]^(1.0 - alpha)

    # get average income and average working hours
    INC[it] = w[it] * LL[it] / workpop[it]
    HH[it] = HH[it] / workpop[it]

    # get difference on goods market
    DIFF[it] = YY[it] - CC[it] - II[it] - GG[it]

end


# function for calculating government parameters
function government(it)

    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # set government quantities and pension payments
    GG[it] = gy * YY[0]
    BB[it] = by * YY[0]
    PP[it] = 0.0

    for ij in JR:JJ
        for ik in 1:SS
            PP[it] = PP[it] + pen[ij, ik, it] * m[ij, ik, it]
        end
    end

    # calculate government expenditure
    expend = GG[it] + (1.0 + r[it]) * BB[it] - (1.0 + n_p) * BB[itp]

    # get budget balancing tax rate
    if (tax[it] == 1)
        tauc[it] = (expend + -(tauw[it] * w[it] * LL[it] + taur[it] * r[it] * AA[it])) / CC[it]
        p[it] = 1.0 + tauc[it]
    elseif (tax[it] == 2)
        tauw[it] = (expend + -tauc[it] * CC[it]) / (w[it] * LL[it] + r[it] * AA[it])
        taur[it] = tauw[it]
    elseif (tax[it] == 3)
        tauw[it] = (expend + -(tauc[it] * CC[it] + taur[it] * r[it] * AA[it])) / (w[it] * LL[it])
    else
        taur[it] = (expend + -(tauc[it] * CC[it] + tauw[it] * w[it] * LL[it])) / (r[it] * AA[it])
    end

    taxrev[1, it] = tauc[it] * CC[it]
    taxrev[2, it] = tauw[it] * w[it] * LL[it]
    taxrev[3, it] = taur[it] * r[it] * AA[it]
    taxrev[4, it] = sum(taxrev[1:3, it])

    # get budget balancing social security replacement rate
    EP_total = PP[it] / kappa[it]
    taup[it] = kappa[it] * EP_total / (w[it] * LL[it])

end


# function for calculating implicit taxes in the pension system
function implicit_taxes(it)


    for ij in 1:JR-1
        itp = year(it, ij, JJ)
        tau_impl[ij, it] = 0.0
        for ijp in JJ:-1:ij+1
            itp = year(it, ij, ijp)
            if (ijp >= JR)
                tau_impl[ij, it] = tau_impl[ij, it] + kappa[itp] * INC[itp]
            end
            tau_impl[ij, it] = tau_impl[ij, it] / (1.0 + rn[itp])
        end
        tau_impl[ij, it] = taup[it] - (1.0 - lambda[it]) * tau_impl[ij, it] / (Float64(JR - 1) * INC[it])
    end
end



# initializes transitional variables
function initialize_trn()

    println("TRANSITION PATH")

    println("ITER       H     K/Y     C/Y     I/Y       r       w        DIFF")

    #= set up population structure
    for it in 1:TT
        pop[1, it] = (1.0 + n_p) * pop[1, it-1]
        for ij in 2:JJ
            pop[ij, it] = pop[ij-1, it-1]
        end
    end

    for it in 1:TT
        for ij in 1:JJ
            m[ij, it] = pop[ij, it] / pop[1, it]
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
        rn[it] = r[it] * (1.0 - taur[it])
        w[it] = w[0]
        wn[it] = w[it] * (1.0 - tauw[it] - taup[it])
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

        pen[:, :, it] = pen[:, :, 0]
        penp[:, :, it, :] = penp[:, :, 0, :]
        taup[it] = taup[0]
        PP[it] = PP[0]
        taxrev[:, it] = taxrev[:, 0]

        c_coh[:, :, it] = c_coh[:, :, 0]
        l_coh[:, :, it] = l_coh[:, :, 0]
        y_coh[:, :, it] = y_coh[:, :, 0]
        a_coh[:, :, it] = a_coh[:, :, 0]

        aplus[:, :, :, :, :, :, it] = aplus[:, :, :, :, :, :, 0]
        epplus[:, :, :, :, :, :, it] = epplus[:, :, :, :, :, :, 0]
        c[:, :, :, :, :, :, it] = c[:, :, :, :, :, :, 0]
        l[:, :, :, :, :, :, it] = l[:, :, :, :, :, :, 0]
        phi[:, :, :, :, :, :, it] = phi[:, :, :, :, :, :, 0]
        VV[:, :, :, :, :, :, it] = VV[:, :, :, :, :, :, 0]
        RHS[:, :, :, :, :, :, it] = RHS[:, :, :, :, :, :, 0]
    end

end



# computes the transition path of the economy
function get_transition()


    # initialize remaining variables
    if (!lsra_on)
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

        # get implicit tax rates
        for it in 1:TT
            implicit_taxes(it)
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
        if lsra_on
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
            if (abs(DIFF[it] / YY[it]) * 100.0 < sig)
                itcheck = itcheck + 1
            end
        end

        #itmax = maxloc(abs(DIFF[1:TT]/YY[1:TT]), 1)
        itmax = argmax(abs.(DIFF[1:TT] ./ YY[1:TT]))

        # check for convergence and write screen output
        if (!lsra_on)

            check = iter > 1 && itcheck == TT-1 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0

            fmt = "%4i " * "%8.2f "^6 * "%12.5f\n" #'(i4,6f8.2,f12.5)'

            @eval @printf $fmt $iter HH[TT] ([5.0*KK[TT], CC[TT], II[TT]]/YY[TT]*100)... ((1.0+r[TT])^0.2-1.0)*100.0 w[TT] DIFF[TT]/YY[TT]*100.0

        else
            check = iter > 1 && itcheck == TT-1 && lsra_comp/lsra_all > 0.99999 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0

            fmt = "%4i " * "%12.5f"^3 * "\n" # '(i4,3f12.5)'

            @eval @printf $fmt $iter  (lsra_comp/lsra_all*100.0) (Lstar^(1.0/(1.0-1.0/gamma))-1.0)*100.0 DIFF[$itmax]/YY[$itmax]*100.0
        end

        # check for convergence
        if check
            #call toc
            for it in 1:TT
                if (!lsra_on)
                    output(it)
                end
            end
            output_summary()
            return
        end
    end

    #call toc
    for it in 1:TT
        if (!lsra_on)
            output(it)
        end
    end
    output_summary()
    println("ATTENTION: NO CONVERGENCE ###")

    #write(*,'(a/)')'ATTENTION: NO CONVERGENCE ###'

end


# function for writing output
function output(it)

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
        for ik in 1:SS
            for ia = 0:NA
                for ir = 0:NR
                    for ip = 1:NP
                        for is = 1:NS
    
                            # consumption
                            if (c[ij, ik, ia, ir, ip, is, it] > 0.0)
                                temp = log(c[ij, ik, ia, ir, ip, is, it])
                                exp_c[ij] = exp_c[ij] + temp*phi[ij, ik, ia, ir, ip, is, it]
                                var_c[ij] = var_c[ij] + temp^2*phi[ij, ik, ia, ir, ip, is, it]
                                mas_c[ij] = mas_c[ij] + phi[ij, ik, ia, ir, ip, is, it]
                            end
    
                            if (l[ij, ik, ia, ir, ip, is, it] > 0.01)
    
                                # hours
                                temp = log(l[ij, ik, ia, ir, ip, is, it])
                                exp_l[ij] = exp_l[ij] + temp*phi[ij, ik, ia, ir, ip, is, it]
                                var_l[ij] = var_l[ij] + temp^2*phi[ij, ik, ia, ir, ip, is, it]
                                mas_l[ij] = mas_l[ij] + phi[ij, ik, ia, ir, ip, is, it]
    
                                # earnings
                                temp = log(w[it]*eff[ij, ik]*theta[ip]*eta[is]*l[ij, ik, ia, ir, ip, is, it])
                                exp_y[ij] = exp_y[ij] + temp*phi[ij, ik, ia, ir, ip, is, it]
                                var_y[ij] = var_y[ij] + temp^2*phi[ij, ik, ia, ir, ip, is, it]
                                mas_y[ij] = mas_y[ij] + phi[ij, ik, ia, ir, ip, is, it]
                            end
                        end
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
     @eval @printf file_output $fmt_val taup[$it]*w[$it]*LL[$it] mean(pen[JR,:, 0]) PP[$it]
     @eval @printf file_output $fmt_rate "(in %) " ([taup[$it], kappa[$it], PP[$it]/YY[$it]]*100.0)... 
 
     fmt_val = " "^8*"%8.2f "^2*"\n"
     fmt_rate = "%s "*"%8.2f "^2*"\n\n"
     @printf file_output "%s \n" "LSRA          SV      BA"
     @eval @printf file_output $fmt_val SV[$it] BA[$it] 
     @eval @printf file_output $fmt_rate "(in %) " ([SV[$it], BA[$it]]/YY[$it]*100.0)...
 
 
     # check for the maximium grid point used
     iamax, iemax = check_grid(it)
 
     itp = year(it, 1, 2)
 
     @printf file_output "%s\n" "   IJ      CONS     LABOR  EARNINGS    INCOME    INCTAX      PENS    ASSETS    VAR(C)    VAR(L)    VAR(Y)      LSRA     VALUE     FLC       IAMAX     IEMAX"
 
     fmt = "%3i "*"%10.3f "^13*"%10i %10i \n"
 
     for ij in 1:JJ
         @eval @printf file_output $fmt $ij  sum(c_coh[$ij,:,$it])/INC[0]  sum(l_coh[$ij,:, $it]) ([w[$it]*sum(y_coh[$ij, :, $it]), wn[$it]*sum(y_coh[$ij, :, $it])+rn[$it]*sum(a_coh[$ij, :, $it]), tauw[$it]*w[$it]*sum(y_coh[$ij, :, $it])+taur[$it]*r[$it]*sum(a_coh[$ij, :, $it]), sum(pen[$ij, :, $it]-taup[$it]*w[$it]*y_coh[$ij, :, $it]) , 5.0*sum(a_coh[$ij, :, $it])]/INC[0])... $var_c[$ij] $var_l[$ij] $var_y[$ij] sum(v_coh[$ij,:, $it]) sum(VV_coh[$ij,:,$it]) sum(FLC[$ij, :, $it]) $iamax[$ij] $iemax[$ij]
     end
 
     @printf file_output "%s \n\n" "--------------------------------------------------------------------"
 
 end 



# subroutine that checks for the maximum gridpoint used
function check_grid(it)

    iamax = zeros(JJ)
    iemax = zeros(JJ)

    for ij = 1:JJ
        for ik = 1:SS
            # check for the maximum asset grid point used at a certain age
            for ia = 0:NA
                for ir = 0:NR
                    for ip = 1:NP
                        for is = 1:NS
                            if (phi[ij, ik, ia, ir, ip, is, it] > 1e-6)
                                iamax[ij] = ia
                            end
                        end
                    end
                end
            end
        

            # check for the maximum earning point grid point used at a certain age
            for ir = 0:NR
                for ia = 0:NA
                    for ip = 1:NP
                        for is = 1:NS
                            if (phi[ij, ik, ia, ir, ip, is, it] > 1e-6)
                                iemax[ij] = ir
                            end
                        end
                    end
                end
            end
        end
    end
    
    return iamax, iemax
end 


# writes summary output
function output_summary()

    HEV = OffsetArray(zeros( length(-(JJ-2):TT) ) , -(JJ-2):TT )
    mas = OffsetArray(zeros( length( -(JJ-2):0)), -(JJ-2):0)

    # aggregate ex post welfare changes of current generations
    HEV .= 0.0
    mas .= 0.0

    for ij = JJ:-1:2
        for ik = 1:SS
            for ia = 0:NA
                for ir = 0:NR
                    for ip = 1:NP
                        for is = 1:NS
                            if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                                continue
                            end
                            HEV_help = ((VV[ij, ik, ia, ir, ip, is, 1]/max(VV[ij, ik, ia, ir, ip, is, 0], -1e10))^(1.0/egam)-1.0)*100.0
                            HEV[-(ij-2)] = HEV[-(ij-2)] + HEV_help*phi[ij, ik, ia, ir, ip, is, 1]
                            mas[-(ij-2)] = mas[-(ij-2)] + phi[ij, ik, ia, ir, ip, is, 1]
                        end
                    end
                end
            end
        end
    end

    HEV[-(JJ-2):0] = HEV[-(JJ-2):0]./parent(mas)

    # calculate ex ante welfare of future generations
    for it = 1:TT
        for ik = 1:SS
            HEV[it] = ((VV_coh[1, ij, it]/VV_coh[1, ik, 0])^(1.0/egam)-1.0)*100.0
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


# function for calculating lsra payments
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
        for ik = 1:SS
            for ia = 0:NA
                for ir = 0:NR
                    for ip = 1:NP
                        for is = 1:NS

                            # for not for anything for an agent at retirement without pension and savings
                            if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10) )
                                v[ij, ik, ia, ir, ip, is, 1] = 0.0
                                continue
                            end

                            # get today's utility
                            VV_1 = VV[ij, ik, ia, ir, ip, is, 1]

                            # get target utility
                            VV_0 = VV[ij, ik, ia, ir, ip, is, 0]


                            # get derivative of the value function
                            dVV_dv = margu(c[ij, ik, ia, ir, ip, is, 1],l[ij, ik, ia, ir, ip, is, 1], 1)


                            # calculate change in transfers
                            v_tilde = (VV_0-VV_1)/dVV_dv

                            # restrict z_tilde to income maximum
                            v_tilde = max(v_tilde, -((1.0+rn[1])*a[ia] + penp[ij, ik, 1, ir] + wn[1]*eff[ij, ik]*theta[ip]*eta[is]*0.99 + v[ij, ik, ia, ir, ip, is, 1]))

                            # check whether individual is already compensated
                            lsra_all = lsra_all + phi[ij, ik, ia, ir, ip, is, 1]*m[ij, ik, 1]

                            if (abs((VV_1-VV_0)/VV_0)*100.0 < sig) 
                                lsra_comp = lsra_comp + phi[ij, ik, ia, ir, ip, is, 1]*m[ij, ik, 1]
                            end

                            # calculate total transfer
                            v[ij, ik, ia, ir, ip, is, 1] = v[ij, ik, ia, ir, ip, is, 1] + damp*v_tilde

                            # aggregate transfers by cohort
                            v_coh[ij, ik, 1] = v_coh[ij, ik, 1] + v[ij, ik, ia, ir, ip, is, 1]*phi[ij, ik, ia, ir, ip, is, 1]

                        end
                    end
                end
            end
        end
    end

    # aggregate transfers in year 1
    for ij = 2:JJ
        for ik = 1:SS
            SV[1] = SV[1] + v_coh[ij, ik, 1]*m[ij, ik, 1]
        end
    end

    # initialize present value variables
    PV_t = 0.0
    PV_0 = 0.0
    PV_trans = 0.0

    # calculate present value of utility changes (in monetary values)
    PV_t_skills = zeros(SS)
    PV_trans_skills = zeros(SS)
    PV_0_skills = zeros(SS)
    
    for it = TT:-1:1
        for ik = 1:SS
            # get today's ex ante utility
            EVV_t = damp*VV_coh[1, ik, it]

            # get damped target utility
            EVV_0 = damp*VV_coh[1, ik, 0]

            # get derivative of expected utility function
            dEVV_dv = 0.0
            for ip = 1:NP
                for is = 1:NS
                    dEVV_dv = dEVV_dv + margu(c[1, ik, 0, 0, ip, is, it], l[1, ik, 0, 0, ip, is, it], it)*phi[1, ik, 0, 0, ip, is, it]
                end
            end

            # calculate present values
            if (it == TT)
                PV_t_skills[ik]     = EVV_t/dEVV_dv    *(1.0+r[it])/(r[it]-n_p)
                PV_0_skills[ik]     = EVV_0/dEVV_dv    *(1.0+r[it])/(r[it]-n_p)
                PV_trans_skills[ik] = v[1, ik, 0, 0, 1, 1, it]*(1.0+r[it])/(r[it]-n_p)
            else
                PV_t_skills[ik]     = PV_t_skills[ik]    *(1.0+n_p)/(1.0+r[it+1]) + EVV_t/dEVV_dv
                PV_0_skills[ik]     = PV_0_skills[ik]    *(1.0+n_p)/(1.0+r[it+1]) + EVV_0/dEVV_dv
                PV_trans_skills[ik] = PV_trans_skills[ik]*(1.0+n_p)/(1.0+r[it+1]) + v[1, ik, 0, 0, 1, 1, it]
            end
        end
    end


    # calculate the constant utility gain/loss for future generations
    Lstar = (sum(PV_t_skills)-sum(PV_trans_skills)-SV[1])/sum(PV_0_skills)

    # calculate compensation payments for future cohorts
    for it = TT:-1:1
        for ik = 1:SS
            # get today's ex ante utility
            EVV_t = damp*VV_coh[1, ik, it]

            # get target utility
            EVV_0 = damp*VV_coh[1, ik, 0]*Lstar

            # get derivative of expected utility function
            dEVV_dv = 0.0
            for ip = 1:NP
                for is = 1:NS
                    dEVV_dv = dEVV_dv + margu(c[1, ik, 0, 0, ip, is, it], l[1, ik, 0, 0, ip, is, it], it)*phi[1, ik, 0, 0, ip, is, it]
                end
            end

            # compute change in transfers (restricted)
            v_tilde = (EVV_0-EVV_t)/dEVV_dv

            # calculate cohort transfer level
            v[1, ik, 0, 0, :, :, it] = v[1, ik, 0, 0, :, :, it] .+ v_tilde

            # aggregate transfers
            v_coh[1, ik, it] = v[1, ik, 0, 0, 1, 1, it]
            SV[it] = SV[it] + v_coh[1, ik, it]*m[1, ik, it]
        end

    end

    # determine sequence of LSRA debt/savings
    BA[2] = SV[1]/(1.0+n_p)
    for it = 3:TT
        BA[it] = ((1.0+r[it-1])*BA[it-1] + SV[it-1])/(1.0+n_p)
    end

end 
