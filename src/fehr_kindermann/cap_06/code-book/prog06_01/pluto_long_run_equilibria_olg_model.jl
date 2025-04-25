### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ b8dd144f-2b43-40e6-af8d-db549a3387ea
begin
	using Printf
	using NLsolve
end

# ╔═╡ 29cebbb7-6986-4d47-96ad-0ccebb993dfc
html"<button onclick='present()'>present</button>"

# ╔═╡ 63da7e09-da7d-441a-b967-892d7dcf1f78
md"""

# Introduction to Computational Economics Using Fortran
"""

# ╔═╡ 75ba6402-0628-11ee-02ad-bbd3ec9278e3
md"""
# Cap 6: The Overlapping generations models

* Individual saving behaviour endogenously determines the capital stock.


"""

# ╔═╡ 723a3b77-07f3-4ca2-839f-fa25f5034198
md"""
# General structure and long-run equilibrium

* Demographics.
* Households' decisions.
* Firms' decisions.
* The goverment.
* Market equilibrium
"""

# ╔═╡ 33d59c93-b8bb-4f37-b56a-b3dc9ba463f0
md"""
# Demographics
* Households in the model live for three periods.
* In each succesive period $t$ a new cohort is born, where $N_t$ in this cohort grows at a rate $n_{p,t}$

```math
	N_t = (1+n_{p,t})N_{t-1}
```
* In period $t$ the total population size is $N_t + N_{t-1} + N_{t-2}$
* Assume that households work in the first two periods of their life and retire in the third, $N_t + N_{t-1}$ defines the size of the workforce in $t$.
"""

# ╔═╡ fde05f36-4005-449b-820a-54c13efaff45
md"""
# Households' decisions.

* They have to pay taxes on income and consumption as well as a payroll tax to the pension system, in reward for which ther receive pension benefits $pen_{t+2}$ in the third period of life.
* The **dynamic budget constraints** in the three diferent periods are given by

```math
	p_{t}c_{1,t} = w_{t}^{n} - a_{2,t+1}
```
```math
	p_{t+1}c_{2,t+1} = R_{t+1}^{n} a_{2,t+1} +w_{t+1}^{n} - a_{3,t+2}
```

```math
	p_{t+2}c_{3,t+2} = R_{t+2}^{n} a_{3,t+2} + pen_{t+2}
```

Where

```math
	w_t^n = w_t(1-\tau_{t}^{w}-\tau_{t}^{p})
```

```math
	R_{t}^{n} = 1+ r_t(1-\tau_t^r)
```

```math
	p_t = 1 + \tau_r^c
```

define net wage, net interest rate, and consumer prices, respectively. $\tau_{t}^{c}, \tau_{t}^{w}, \tau_{t}^{p}, \tau_{t}^{r}$ denote consumption, labor-income, payroll, and capital-income tax rates.
"""

# ╔═╡ 255f6391-0241-44a9-b89c-acc223ea57df
md"""
# Utility functions

* Households born in period $t$ maximize the utility function

```math
	U_{t} = U(C_{1,t}, C_{2,t}, C_{3,t}) = u(c_{1,t}) + \beta u(c_{2,t+1}) + \beta^2 u(c_{3,t+2})
```

subject to the **intertemporal budget constraint**

```math
	p_{t}c_{1,t} + \dfrac{p_{t+1}c_{2,t+1}}{R_{t+1}^n} +  \dfrac{p_{t+2}c_{3,t+2}}{R_{t+1}^n R_{t+2}^n} = w_t^n + \dfrac{w_{t+1}^n}{R_{t+1}^n} + \dfrac{pen_{t+2}}{R_{t+1}^nR_{t+2}^n} =: W_{1,t}
```

"""

# ╔═╡ a9935a8b-1f29-41c1-b8a9-8b5e3821ba3b
md"""
# Agreggation

* Since different households are living in the same period $t$, the demand for goods and the supply of labour and capital need to be determined by aggregrating individual households decisions.

* Normalize all aggregate variables to the size of the newborn cohort.

* Size each working households is assumed to supply one unit of labour inelastically, aggregate labour supply in period $t$ can be computed from

```math
	L_t = \dfrac{1}{N_t} (N_t + N_{t-1}) = 1 + \dfrac{1}{1+n_{p,t}} = \dfrac{2+n_{p,t}}{1+n_{p,t}}
```

* We derive aggregate per capita assets
```math
	A_t = \dfrac{1}{N_t} (N_{t-1}a_{2,t} + N_{t-2}a_{3,t})  = \dfrac{a_{2,t}}{1+n_{p,t}} + \dfrac{a_{3,t}}{(1+n_{p,t})(1+n_{p,t-1}) }
```

where $a_{2,t}$ and $a_{3,t}$ denote savings of the two older cohorts built in the previous year, and per capita consumption from

```math
	C_t = \dfrac{1}{N_t}(N_tC_{1,t} + N_{t-1}C_{2,t} + N_{t-2}C_{3,t})) = c_{1,t} +  \dfrac{c_{2,t}}{1+n_{p,t}} + \dfrac{c_{3,t}}{(1+n_{p,t})(1+n_{p,t-1}) }
```

"""

# ╔═╡ b1ead0ce-c090-4043-8584-37196333a991
md"""
# Firms' decisions
* Production technology
```math
Y_t = f(K_t, L_t) \qquad \text{or} \qquad y_t = \dfrac{K_t}{L_t} = f(k_t)
```

where $y_t$ denote per capita output and $k_t = \dfrac{K_t}{L_t}$ capital intensity.

* Firms are assumed to maximize (per capita) profits $\pi_t = y_t - (r_t+\delta)k_t-w_t$ by choosing the optimal capital intensity.

* Firms produce under perfect competition conditions. Therefore equilibrium factor prices $r_t$ and $w_t$ are derived from FOC and zero-profit condition as

```math
	r_t = f'(k_t) - \delta
```
```math
	w_t = f(k_t) - f'(k_t)k_t
```
"""

# ╔═╡ 724d2ab8-a59d-4774-ae0e-5d6a88f33dbc
md"""
# The goverment

* The goverment sector comprisses two separate budges.
* It finances the provision of a public good $G_t$ and interest payment on public debt $r_tB_t$ by mean of tax reveneus and deficit spending

```math
	T_t + (1+n_{p,t+1})B_{t+1} - B_{t} = G_t + r_tB_t
```

where 
```math
	T_t = \tau_{t}^{c}C_t + \tau_{t}^{w} w_tL_t + \tau_{t}^{r}r_tA_t
```

define tax reveneus from consumption and income taxation. The amount of public good is computed depending on the population structure

```math
	G_t = \dfrac{1}{N_t}(N_tg_1 + N_{t-1} g_2 + N_{t-2}g_3)
```

```math
	 = g_1 +  \dfrac{g_{2}}{1+n_{p,t}} + \dfrac{g_{3}}{(1+n_{p,t})(1+n_{p,t-1}) }
```

where $g_j$ represent per capita coeficient of cohort $j$ for public good provision. Furthemore, we determine the deficit path exogenously such that

```math
	 B_{t+1} = b_{y,t} Y_{t+1}
```

holds in every period. 

Given expenditure for the public good $G_t$ and the deficit path $b_{y,t}$, we let of the three tax rates balance the periodical budget of the goverment.
"""

# ╔═╡ 26359fa6-c694-4cd7-afb3-d275fa3a30cf
md"""
# The goverment

With respect to the second goverment budget we assume that in each period $t$ the level of pension benefit is computed as a fraction $\kappa_t$ of the previous gross wage

```math
	 pen_{t} = \kappa_t w_{t-1}
```

The pension system is fully pay-as-you-go financed. Hence, contributions from workers have to finance benefit payments to retirees. 

```math
	 \tau_t^p w_t(N_t+N_{t-1}) = pen_t N_{t-2} \quad \text{or} \quad \tau_{t}^pw_t=\dfrac{pen_t}{(2+n_{p,t})(1+n_{p,t-1})}
```

"""

# ╔═╡ cc88ac8e-3da3-48bd-8d39-efdc4750d310
md"""
# Market equilibrium

For a market equilibrium the goods market and factor markets need to be balanced. 

* On the capital market, of a closed economy, interest rates are determined so as to equalize households' aggregate saving $A_t$, with the demand for capital of the firms $K_t$ and the goverment $B_t$, *i.e.*


```math
	 A_t = K_t + B_t
```

* In the labor market, wage are set in a way that guarantees that aggregate labour supply of households match labour demand of the firms.
* Finally, in the goods market total output of the firms must be used either for private and public consumption or for holding the new capital stock, *i.e.*

```math
	 Y_t = C_t + G_t + (1 + n_{p,t+1})K_{t+1} - (1-\delta)K_t
```

"""

# ╔═╡ 3262b340-4059-496d-986d-7aa6983c924e
md"""
# Computation of the long-run equilibrium

Steady states are characterized by the fact that per capita variables are constant over time.

## Functional forms

* We have to define concrete functional forms of preferences and production functions.
* The CRRA utility function given by

```math
	 u(c) = \dfrac{c^{1-\frac{1}{\gamma}}}{1-\frac{1}{\gamma}} \quad \text{with} \quad u'(c) = c^{-\dfrac{1}{\gamma}}
```

where $\gamma$ denotes the intertemporal elasticity of sustitution between consumption in different years.

* The production technology is assumed to be a Cobb-Douglas type
```math
	 f(k) = K^\alpha 
```

* Interest and wage rates are then determined by
```math
	 r = \alpha K^{\alpha-1}-\delta \qquad \text{and} \qquad w = (1-\alpha)K^\alpha 
```

"""

# ╔═╡ 55b3248c-6990-43f9-8d8f-234a7fe6d1ab
md"""

# Exogenous parameters

* Endogenous
  * Households decisions
  * Aggregate quantities
  * Factor prices
  * Payroll tax rate
  * One aditional tax rate (consumption, labour income, or capital income which balances the goverment budget)

* Exogenous
  * Preferences parameters ($\gamma, \beta$).
  * Technology parameters ($\alpha, \delta$).
  * Population growth rate, $n_p$
  * Goverment parameters ($g_j, \tau^w, \tau^r, b_y$).
  * Replacement rate $\kappa$ of the pension system.

"""

# ╔═╡ 63158864-d297-4064-a8a9-40c959f009fb
md"""
# Code
"""

# ╔═╡ f3ac39d3-5732-49ce-8d00-9b4caaa5dafb
begin	
	# model parameters
	global γ = 0.5
	global egam = 1-(1/γ)
	global β = 0.9
	global α = 0.3
	global δ = 0.0
	global by = 0
	global κ = 0.0
	global nₚ = 0.2
	global tax = 1
	
	global g = [0.12, 0.12, 0.0]
	
	# solve the steady state equation system
	
	function eqns(F, x)
	    # model variables
	    global w = 0.0
	    global r = 0.0
	    global wn = 0.0
	    global Rₙ = 0.0
	    global p = 0.0
	    global τw = 0.0
	    global τr = 0.0
	    global τc = 0.0
	    global τp = 0.0
	    global pen = 0.0
	
	    global K = 0.0
	    global L = 0.0 
	    global Y = 0.0
	    global A = 0.0
	    global C = 0.0
	    global B = 0.0
	    global G = 0.0
	    global I = 0.0
	
	    global a = [0.0, 0.0, 0.0]
	    global c = [0.0, 0.0, 0.0] 
	
	    #initialize labor supply, pension payments and tax rates
	
	    global L = (2.0+nₚ)/(1.0+nₚ)
	    global τp = κ/((2.0+nₚ)*(1.0+nₚ))
	    global τc = 0.0
	    global τw = 0.0
	    global τr = 0.0
	
	    K = x[1]
	    
	    if tax == 1
	        τc = x[2]
	    elseif tax==2
	        τw = x[2]
	        τr = x[2]
	    elseif tax == 3
	        τw = x[2]
	    else tax > 3
	        τr = x[2]
	    end
	    
	    # factor prices and pension payments
	    r = α*(K/L)^(α-1.0)-δ
	    w = (1.0-α)*(K/L)^α
	    wn = w*(1.0-τw-τp)
	    Rₙ = 1.0+r*(1.0-τr)
	    p = 1.0+τc
	    pen = κ*w
	
	    # individual decisions
	    PVI = wn + wn/Rₙ + pen/Rₙ^2
	    PSI = p*(1.0 + (β*Rₙ)^γ/Rₙ + (β*Rₙ)^(2*γ)/Rₙ^2)
	
	    c[1] = PVI/PSI
	    c[2] = (β*Rₙ)^γ*c[1]
	    c[3] = (β*Rₙ)^γ*c[2]
	    a[1] = 0.0
	    a[2] = wn - p*c[1]
	    a[3] = wn + Rₙ*a[2] - p*c[2]
	
	    # quantities
	    Y = K^α * L^(1.0-α)
	    C = c[1] + c[2]/(1.0+nₚ) + c[3]/(1.0+nₚ)^2
	    G = g[1] + g[2]/(1.0+nₚ) + g[3]/(1.0+nₚ)^2
	    A = a[2]/(1.0+nₚ) + a[3]/(1.0+nₚ)^2
	    I = (nₚ+δ)*K
	    B = by*Y
	
	    # get equations defining general equilibrium
	    F[1] = K + B - A
	    F[2] = τc*C + τw*w*L + τr*r*A - (r-nₚ)*B - G
	end
		
	x = nlsolve(eqns, [ 0.1; 0.1])


	function calculate_household_utility()
		# Calculate household utility
		global util = 0
		
		for i in 1:3
		    util +=  β^(i-1)*c[i]^egam/egam
		end
		return util
	end

	util = calculate_household_utility()
end

# ╔═╡ 6dc76e70-fae6-4f73-9e38-afc3feb6c330

begin

	@printf("Steady state equilibrium\n")
	@printf("     c1     c2     c3      Y      w      r      U\n")
	@printf("   %0.2f   %0.2f   %0.2f    %0.2f   %0.2f   %0.2f    %0.2f   \n\n\n", c[1], c[2], c[3], Y, w, r,util)


	@printf("     a2     a3      K \n")
	@printf("   %0.2f   %0.2f   %0.2f\n\n\n", a[2], a[3], K)
 
	@printf("      Y      C      G      I      DIFF\n")
	@printf("   %0.2f   %0.2f   %0.2f    %0.2f   %0.2f \n\n\n", Y, C, G, I, Y - C - G - I)


end


# ╔═╡ 888aa006-6d0a-48dd-93d8-7d5ba0c7abf3
md"""
## Results
|    |     by |   kappa |   np |   tauw |   tawr |   tauc |   c1 |   c2 |   c3 |    K |      U |
|---:|-------:|--------:|-----:|-------:|-------:|-------:|-----:|-----:|-----:|-----:|-------:|
|  0 |  0     |     0   |  0.2 |   0    |   0    |   0.29 | 0.22 | 0.3  | 0.42 | 0.27 |  -9.54 |
|  1 |  0     |     0   |  0.2 |   0.24 |   0.24 |   0    | 0.19 | 0.26 | 0.37 | 0.18 | -10.94 |
|  2 |  0     |     0   |  0.2 |   0    |   0.71 |   0    | 0.28 | 0.3  | 0.33 | 0.27 |  -9.03 |
|  3 |  0     |     0   |  0.2 |   0.37 |   0    |   0    | 0.15 | 0.23 | 0.37 | 0.14 | -12.93 |
|  4 | -0.059 |     0   |  0.2 |   0.22 |   0    |   0    | 0.22 | 0.3  | 0.42 | 0.27 |  -9.54 |
|  5 |  0     |     0.5 |  0.2 |   0    |   0    |   0.37 | 0.14 | 0.23 | 0.37 | 0.14 | -13.01 |
|  6 |  0.099 |     0   |  0.2 |   0    |   0    |   0.6  | 0.14 | 0.23 | 0.37 | 0.14 | -13.04 |
|  7 |  0     |     0   |  0   |   0    |   0    |   0.24 | 0.25 | 0.32 | 0.43 | 0.4  |  -8.74 |

"""

# ╔═╡ 67fdcee5-91ee-42d5-bcb7-6b944c134f31
md"""
# Extensions

* Accounting for variable labour supply.
* Human capital (Education investment).
* Introducing private annuity markets.
  * Uncertain lifespan.
  * Rising life expectacy.
* Cap 11. Dynamic macro II: The stochastic OLG model.
"""

# ╔═╡ 82933d61-8af8-496c-95a3-36288bdbfa24


# ╔═╡ 6bbf0754-994f-4855-ab5d-a9396b646fd4


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
NLsolve = "~4.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "36abd45079f89b81d0b77b20980be50a7ec4bc1e"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d3f758863a47ceef2248d136657cb9c033603641"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.8"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─29cebbb7-6986-4d47-96ad-0ccebb993dfc
# ╟─63da7e09-da7d-441a-b967-892d7dcf1f78
# ╟─75ba6402-0628-11ee-02ad-bbd3ec9278e3
# ╟─723a3b77-07f3-4ca2-839f-fa25f5034198
# ╟─33d59c93-b8bb-4f37-b56a-b3dc9ba463f0
# ╟─fde05f36-4005-449b-820a-54c13efaff45
# ╟─255f6391-0241-44a9-b89c-acc223ea57df
# ╟─a9935a8b-1f29-41c1-b8a9-8b5e3821ba3b
# ╟─b1ead0ce-c090-4043-8584-37196333a991
# ╟─724d2ab8-a59d-4774-ae0e-5d6a88f33dbc
# ╟─26359fa6-c694-4cd7-afb3-d275fa3a30cf
# ╟─cc88ac8e-3da3-48bd-8d39-efdc4750d310
# ╟─3262b340-4059-496d-986d-7aa6983c924e
# ╟─55b3248c-6990-43f9-8d8f-234a7fe6d1ab
# ╟─63158864-d297-4064-a8a9-40c959f009fb
# ╠═b8dd144f-2b43-40e6-af8d-db549a3387ea
# ╠═f3ac39d3-5732-49ce-8d00-9b4caaa5dafb
# ╠═6dc76e70-fae6-4f73-9e38-afc3feb6c330
# ╟─888aa006-6d0a-48dd-93d8-7d5ba0c7abf3
# ╟─67fdcee5-91ee-42d5-bcb7-6b944c134f31
# ╠═82933d61-8af8-496c-95a3-36288bdbfa24
# ╠═6bbf0754-994f-4855-ab5d-a9396b646fd4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
