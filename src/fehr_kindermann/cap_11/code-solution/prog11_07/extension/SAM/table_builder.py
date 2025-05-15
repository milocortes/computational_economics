import pypst
import pandas as pd
import numpy as np 

sam = pd.read_csv("outputs/sam_olg.csv")

sam = sam.set_index("Unnamed: 0")
sam.index.names = ['']

sam["Total"] = sam.replace(np.nan, 0.0).sum(axis = 1)
sam.loc["Total"] = sam.replace(np.nan, 0.0).sum(axis = 0)
sam["Balance"] = sam["Total"] - sam.loc["Total"]

table = pypst.Table.from_dataframe(sam.round(3).replace({np.nan : ""}))
#table.stroke = "(x: none)"
table.align = "center + horizon"

figure = pypst.Figure(table, caption='"SAM. Values from the Modelo"')

with open("typst/sam_values.typ", mode="wt") as f:
    f.write(figure.render())

