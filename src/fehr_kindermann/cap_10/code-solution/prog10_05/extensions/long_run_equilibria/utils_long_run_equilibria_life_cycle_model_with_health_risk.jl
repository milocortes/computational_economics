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

