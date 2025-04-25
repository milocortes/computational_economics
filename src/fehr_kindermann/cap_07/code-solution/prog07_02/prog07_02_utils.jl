# initializes variables and government parameters
function initialize()

    # set model parameters
    n_p[0] = 0.2
    n_p[1:TT] .= 0.2
    gy = 0.195
#   gy = 0d0
    lsra_on = false
    smopec = false

    # set reform values

    @. begin
        by[0:TT] = 0.0
        kappa[1:TT] = 0.5
        eps[1:TT] = 0
        tauk[1:TT] = 0.0
        tax[1:TT] = 1
    end

    tax[0] = 1
    tauk[0] = 0.0
    kappa[0] = 0.0
    eps[0] = 0

    # initialize tax rates shadow wages and pensions
    @. begin
        tauc = 0.0
        tauw = 0.0
        taur = 0.0
        taup = 0.0
        pen = 0.0
        mu = 0.0
    end

    if  (nu > 0)
        mu[JJ, :, :] .= 0.5
    end

    # initialize assets, LSRA payments and debt holdings

    @.begin
        a = 0.0
        c = 0.0
        l = 0.0
        v = 0.0
        YY = 0.0
        BA = 0.0
        BF = 0.0
        TB = 0.0
        TXR = 0.0
    end
    # human capital profile
    @. begin
        h[:, 1] = 1.0
        h[:, 2] = 2.0
        h[JR:JJ, :] = 0.0
    end
    # size of cohorts in specific year
    for it in 0:TT
        m[1, 1, it] = 0.2
        m[1, 2, it] = 1.0-m[1, 1, it]
        itm = year(it, 2, 1)
        for ij in 2:JJ
            for ik in 1:SS
                m[ij,ik,it] = m[ij-1,ik,itm]/(1+n_p[it])
            end
        end
    end
end 

# calculates year at which age ij agent is ijp
function year(it, ij, ijp)

    year_val = it + ijp - ij

    if (it == 0 || year_val <= 0)
        year_val = 0
    end
    if (it == TT || year_val >= TT)
        year_val = TT
    end

    return year_val
end 