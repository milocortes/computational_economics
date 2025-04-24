


###############################################################################
# function linint_Gen
#
# For linear interpolation on irregular grids.
###############################################################################
function linint_Gen(x, xi, yi, istart_in)
    ###### ROUTINE CODE #######################################################
    n = length(xi)

    istart = min(max(istart_in, 1), n)

    # if grid value too large, search for first smaller point
    if (xi[istart] > x)
        ial = istart
        while true
            ial = ial - 1
            if (ial <= 1)
                break
            end
            if (xi[ial] <= x)
                break 
            end
        end
        ial = max(ial, 1)
        ial = min(ial, n-1)            
        iar = ial+1
    else
        # if grid value too small, search for first larger point
        iar = istart
        while true 
            iar = iar + 1
            if (iar >= n)
                break
            end
            if (xi[iar] >= x)
                break 
            end
        end
        iar = max(iar, 1)
        iar = min(iar, n)            
        ial = iar-1
    end

    # linearly interpolate between the two points
    phi = 1.0 - (x-xi[ial])/(xi[iar]-xi[ial])        
    return  phi*yi[ial] + (1.0-phi)*yi[iar]

end