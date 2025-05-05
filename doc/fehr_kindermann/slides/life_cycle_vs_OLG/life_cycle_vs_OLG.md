---
title: Introducción a Git
author: Hermilo

theme:
  # Specify it by name for built-in themes
  name: catppuccin-latte
  override:
    default:
      margin:
        percent: 1


---


First Orden Condition
---
<!-- column_layout: [15, 15] -->

<!-- column: 0 -->

The baseline life cycle model
```julia +line_numbers
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

```

<!-- column: 1 -->
 OLG model with survival probabilities
```julia +line_numbers
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
        lab_com = min( max( nu + (1.0-nu)*(a_plus-available)/wage, 0.0) , 1.0-1e-10)
    else
        lab_com = 0.0
    end

    # calculate consumption
    cons_com = max( (available + wage*lab_com - a_plus)/p[it_com] , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = max(varphi*RHS[ij_com+1, ial, ip_com, is_com, itp] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com, itp], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, it_com)^(-gamma) - tomorrow

end 
```
<!-- end_slide -->


Margu
---

<!-- column_layout: [15, 15] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers
# calculates marginal utility of consumption
function margu(cons, lab, im)

    return nu*(cons^nu*(1.0-lab-phi_l*Float64(im))^(1.0-nu))^egam/cons

end 

```


<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers

# calculates marginal utility of consumption
function margu(cons, lab, it)

    return nu*(cons^nu*(1.0-lab)^(1.0-nu))^egam/(p[it]*cons)

end 

```
<!-- end_slide -->

Value Function
---
<!-- column_layout: [15, 15] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers
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
```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers
function valuefunc(a_plus, cons, lab, ij, ip, is, it)

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
        valuefunc = max(varphi*EV[ij+1, ial, ip, is, itp] + (1.0-varphi)*EV[ij+1, iar, ip, is, itp], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^egam/egam + beta*psi[ij+1, itp]*valuefunc

end 

```
<!-- end_slide -->

Solve Households
---
<!-- column_layout: [15, 15] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers
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

```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers
function solve_household(ij_in, it_in)
    global cons_com
    global lab_com

    # get decision in the last period of life
    it = year(it_in, ij_in, JJ)

    #bequest for the olds
    beq[JJ, it] = damp*GAM[JJ, it]*BQ[it] + (1.0-damp)*beq[JJ, it]#

    for ia in 0:NA
        aplus[JJ, ia, :, :, it] .= 0.0
        c[JJ, ia, :, :, it] .= ((1.0+rn[it])*a[ia]+ beq[JJ, it] .+ pen[JJ, it] .+ v[JJ, ia, :, :, it])/p[it]#
        l[JJ, ia, :, :, it] .= 0.0
        VV[JJ, ia, :, :, it] .= valuefunc(0.0, c[JJ, ia, 1, 1, it],l[JJ, ia, 1, 1, it], JJ, 1, 1, it)
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

```
<!-- end_slide -->



Solve Households
---
<!-- column_layout: [15, 15] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers

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
      end
```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers

        for ia in 0:NA

            # determine decision for zero assets at retirement without pension
            if (ij >= JR && ia == 0 && kappa[it] <= 1e-10)
                aplus[ij, ia, :, :, it] = 0.0
                c[ij, ia, :, :, it] = 0.0
                l[ij, ia, :, :, it] = 0.0
                VV[ij, ia, :, :, it] = valuefunc(0.0, 0.0, 0.0, ij, 1, 1, it)
                continue
            end

            for ip in 1:ip_max
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
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
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
                    aplus[ij, ia, :, :, it] .= aplus[ij, ia, 1, 1, it]
                    c[ij, ia, :, :, it] .= c[ij, ia, 1, 1, it]
                    l[ij, ia, :, :, it] .= l[ij, ia, 1, 1, it]
                    VV[ij, ia, :, :, it] .= VV[ij, ia, 1, 1, it]
                end
            end
        end
        # interpolate individual RHS
        interpolate(ij, it)
    end
```
<!-- end_slide -->

Interpolate
---
<!-- column_layout: [10, 10] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers

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

```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers
# for calculating the rhs of the first order condition at age ij
function interpolate(ij, it)

    for ia in 0:NA

        for ip in 1:NP
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is, it] = 0.0
                EV[ij, ia, ip, is, it] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ia, ip, is_p, it],1e-10)
                    lhelp = max(l[ij, ia, ip, is_p, it],1e-10)
                    RHS[ij, ia, ip, is, it] = RHS[ij, ia, ip, is, it] + pi[is, is_p]*margu(chelp, lhelp, it)
                    EV[ij, ia, ip, is, it]  = EV[ij, ia, ip, is, it] + pi[is, is_p]*VV[ij, ia, ip, is_p, it]
                end
                RHS[ij, ia, ip, is, it] = max((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ia, ip, is, it],1e-10)^(-gamma)
                EV[ij, ia, ip, is, it] = (egam*EV[ij, ia, ip, is, it])^(1.0/egam)
            end
        end
    end

end 

```
<!-- end_slide -->

Get Distribution 
---
<!-- column_layout: [10, 10] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers

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
```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers
# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, it] .= 0.0

    # get initial distribution in age 1
    for ip in 1:NP
        phi[1, 0, ip, is_initial, it] = dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is, itm], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi = min(varphi, 1.0)

                    # redistribute households
                    for is_p in 1:NS
                        phi[ij, ial, ip, is_p, it] = phi[ij, ial, ip, is_p, it] + pi[is, is_p]*varphi*phi[ij-1, ia, ip, is, itm]
                        phi[ij, iar, ip, is_p, it] = phi[ij, iar, ip, is_p, it] + pi[is,is_p]*(1.0-varphi)*phi[ij-1,ia,ip,is,itm]
                    end
                end
            end
        end
    end

end 

```
<!-- end_slide -->


Agregation
---
<!-- column_layout: [10, 10] -->

<!-- column: 0 -->
The baseline life cycle model

```julia +line_numbers

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
```

<!-- column: 1 -->
 OLG model with survival probabilities

```julia +line_numbers

# subroutine for calculating quantities in a certain
function aggregation(it)

    # get tomorrow's year
    itp = year(it, 1, 2)
    LL_old = LL[it]
    m_coh = zeros(JJ)

    # calculate cohort aggregates
    c_coh[:, it]  .= 0.0
    l_coh[:, it]  .= 0.0
    y_coh[:, it]  .= 0.0
    a_coh[:, it]  .= 0.0
    VV_coh[:, it] .= 0.0
    m_coh[:]      .= 0.0
    FLC[:,it]     .= 0.0
    beq_coh[:,it] .= 0.0#

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    c_coh[ij, it] = c_coh[ij, it] + c[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    l_coh[ij, it] = l_coh[ij, it] + l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    y_coh[ij, it] = y_coh[ij, it] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    a_coh[ij, it] = a_coh[ij, it] + a[ia]*phi[ij, ia, ip, is, it]
                    # exclude households who dies
                    if (ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        continue
                    end
                    
                    if (aplus[ij, ia, ip, is, it] < 1e-4) 
                        FLC[ij, it] = FLC[ij, it] + phi[ij, ia, ip, is, it]
                    end

                    VV_coh[ij, it] = VV_coh[ij, it] + VV[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    m_coh[ij]      = m_coh[ij] + phi[ij, ia, ip, is, it]
                    beq_coh[ij, it] = beq_coh[ij, it] + a[ia]*(1.0 + rn[it])*(1.0 - psi[ij, it])*phi[ij, ia, ip, is, it]/psi[ij, it] #
                end
            end
        end
    end

```
<!-- end_slide -->