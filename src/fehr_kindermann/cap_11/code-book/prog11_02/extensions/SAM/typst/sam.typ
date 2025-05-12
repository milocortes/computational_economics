
#set page(margin: (
  top: 2cm,
  bottom: 1cm,
  x: 1cm,
))


#align(center,text(24pt)[#smallcaps("The stochastic OLG model with transitional dynamics\nSocial Accounting Matrix")])


#show figure: set block(breakable: true)

#show table.cell.where(y: 0): set text(weight: "bold")


#set page(width: auto, height: auto)

#show table.cell.where(y: 0): strong
#set table(
  stroke: (x, y) => if y == 0 {
    (bottom: 0.7pt + black)
  },
  align: (x, y) => (
    if x > 0 { center }
    else { left }
  )
)

#table(
  columns: 10,
  gutter: 3pt,
  table.header(
    [],
    [*Producción*],
    [*Salarios*],
    [*Ganancias*],
    [*Pensiones*],
    [*Hogares*],
    [*Gobierno (Ingresos)*],
    [*Gobierno (Gasto)*],
    [*Ahorro-Inversión*],
    [*Total*]
  ),
  [Producción], [],[],[],[ ],[$C_t$],[],[$G_t$],[$I_t$],[],
  [Salarios], [$w_t L_t$],[],[],[],[],[],[],[],[],
  [Ganancias], [$(r_t + delta) K_t$],[],[],[],[],[],[],[],[],
  [Pensiones], [],[],[],[],[$w_t L_t tau_t^p$],[],[],[],[],
  [Hogares], [],[$w_t L_t (1-tau_t^w)$],[$r_t A_t (1-tau_t^r)$],[$ overline("pen")_t sum_(j=j_r)^J m_j$ ],[],[],[],[],[],
  [Gobierno (Ingresos)], [],[$w_t L_t tau_t^w$],[$r_t A_t tau_t^r$],[],[$tau_t^c C_t$],[],[],[],[],
  [Gobierno (Gasto)], [],[],[],[],[],[$G_t + (1.0+r_t )B_t – (1.0+n_p)B_(t+1)$],[],[],[],
  [Ahorro-Inversión], [],[],[$r_t A_t$],[],[$r_t A_t tau_t^r$],[],[$(r_t - n_p)B_(t)$],[],[],
  [Total], [],[],[],[],[],[],[],[],[],

)
