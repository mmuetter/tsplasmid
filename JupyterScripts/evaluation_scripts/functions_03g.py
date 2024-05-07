import numpy as np
import pandas as pd

def transition_prop_clear(df, strategy, strain, ab):
    df = df[df.strategy == strategy]
    df = df[df.infected_by_wells.isnull()]
    df = df[df.treatment_with == ab]
    df = df[(df.transferred_phenotype == strain) | (df.turnover_strain == strain) ]
    df["cleared"] = df.phenotype == "UI"
    return np.mean(df["cleared"]), np.std(df["cleared"]), df

def make_contingency_table(df, strain, ab):    
    df = df[df.infected_by_wells.isnull()]
    df = df[df.treatment_with == ab]
    df = df[(df.transferred_phenotype == strain) | (df.turnover_strain == strain) ]
    df["received_strain"] = strain
    df["cleared"] = df.phenotype == "UI"
    cont = df[["strategy", "cleared"]]
    return pd.crosstab(index=cont['strategy'], columns=cont['cleared'])

def load_M(antibiotic):
    ab = "AB"
    path = os.path.join(analysis_path, "obj", "transition_prop", "M_"+antibiotic+".pkl")
    return pd.read_pickle(path).to_numpy()

def create_x_form_dict(dictionary, order):
    x = []
    for i, o in enumerate(order):
        if o in dictionary.keys():
            x.append(dictionary[o])
        else:
            x.append(0)
    return x