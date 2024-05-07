import numpy as np
import pandas as pd

## Count incidence of each Phenotype per Strategy and Replicate
def summerize(Data, groups):
    summary = Data.groupby(groups).count().unstack(fill_value=0).stack()
    summary = summary.reset_index()
    summary = summary.rename(columns={"row": "count", "col" : "fraction"})
    for i, line in summary.iterrows():
        #summary.loc[i, "transfer_n"] = "t" + str(line["transfer_n"]) 
        summary.loc[i, "number_of_patients"] = line["count"]
    return summary