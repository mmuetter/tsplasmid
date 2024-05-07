import pandas as pd
import numpy as np


def dirac(x,a, b):
    dirac = 1/(np.sqrt(np.pi)*a)*np.exp(-abs(x-np.round(x)-b)**2/a**2)
    return dirac

def get_clearence_matrix():
    clearence_rates = pd.read_pickle("obj/Turnover_statistic.pkl")
    clearence_rates = clearence_rates[clearence_rates.time == 1]
    for i, row in clearence_rates.iterrows():
        words = row.history_group.split(" ")
        strain = words[3]
        treatment = " ".join(words[6:])
        clearence_rates.loc[i, "treatment_with"] = treatment
        clearence_rates.loc[i, "phenotype"] = strain 

    clearence_matrix = pd.DataFrame(index=["none", "Mono A", "Mono B", "Combo"], columns = ["UI", "wt", "A_r", "B_r", "AuB", "AB_r"])
    comment = 0
    for strain in clearence_matrix.columns:
        for treatment in clearence_matrix.index:
            if strain == "UI":
                cr = 0
            elif treatment == "none":
                cr = 0
            elif strain == "AB_r":
                cr = 0
            elif strain == "AuB":
                cr = 0
            else:
                cr = clearence_rates[(clearence_rates.phenotype == strain) & (clearence_rates.treatment_with == treatment)]
                cr = cr.fraction.item()
            clearence_matrix.loc[treatment, strain] = cr
    print("There are no values for strain: A\&B and AB. Those values have been replaced by 0. \nComplement this Code in future in Order to implement those rates")
    return clearence_matrix

# The strain is transformed from col-name to index-name
# The strain can only be transformed from a higher to a lower resistance level (otherwise set to 0) (strains that are later more resistant are either unintended infection or conjugation or maybe mutation)
# The strain can not be transformed to UI. This process is already included in clearence (otherwise set to 0)
# Missing values are replaced by 0
# The transformation rates are still very vague
# The groupby.mean method is probably not appropriate to middle over the replicates cause of different numbers of occurence



def get_transformation_rates(antibiotic, comment = 0):
    transformation_rates = pd.read_pickle("obj/stats.pkl")
    transformation_rates = transformation_rates[transformation_rates.treatment_with == antibiotic]

    # The groupby.mean method is probably not appropriate to middle over the replicates cause of different numbers of occurence
    transformation_rates = transformation_rates.groupby(["single_input_strain", "treatment_with", "phenotype"]).mean().reset_index()

    D = pd.DataFrame(index=["UI", "wt", "A_r", "B_r", "A&B", "AB_r"], columns = ["UI", "wt", "A_r", "B_r", "A&B", "AB_r"])
    resistance_level = pd.DataFrame(index=D.index, columns = ["resistance_level"], data=[-1, 0, 1, 1, 2, 3])
    
    ## Norm the fractions to add up to one again
    for in_strain in D.columns:
        norm = transformation_rates[transformation_rates.single_input_strain == in_strain].fraction.sum()
        transformation_rates.loc[transformation_rates.single_input_strain == in_strain, "fraction"] = transformation_rates[transformation_rates.single_input_strain == in_strain].fraction/norm


    for input_strain in D.columns:
        for output_strain in list(D.index):
            if output_strain == "UI":
                D.loc[output_strain, input_strain] = 0
            elif resistance_level.loc[input_strain, "resistance_level"] - resistance_level.loc[output_strain, "resistance_level"] > 0:
                tr = transformation_rates[(transformation_rates.phenotype == output_strain) & (transformation_rates.single_input_strain == input_strain)]
                if tr.empty:
                    D.loc[output_strain, input_strain]  = 0
                    if comment:
                        print("Missing value for input", input_strain, ", output", output_strain, "and antibiotic", antibiotic)
                else:
                    D.loc[output_strain, input_strain] = tr.fraction.item()

            else:
                D.loc[output_strain, input_strain]  = 0
    if comment:
        display(D)
    return D

def get_mu(t, par):
    strategy = par["strategy"]
    start_with_A = par["start_with_A"]
    M = pd.read_excel("Medication.xlsx", index_col=0)
    mu = M.loc[strategy,:]
    if (strategy == "Cycling") & (start_with_A != round(t) % 2 ):
        mu = np.array([0, 1, 0, 0])
    elif (strategy == "Cycling") & (start_with_A == round(t) % 2 ):
        mu = np.array([0, 0, 1, 0])
    return mu

def calc_T(D):
    D = np.matrix(D)
    T = D
    for col in range(6):
        T[col, col] = -D[:, col].sum()
    return T


from random import random
def random_X0():
    X0 = np.array([random(), random(), random(), random(), random(), random()])
    X0 = X0/X0.sum()
    return X0

def reorganize_results_to_dataframe(Xt, strategy, method):
    df = pd.DataFrame(columns = ["UI", "wt", "A_r", "B_r", "A&B", "AB"], data=Xt.y.T)
    df["transfer_n"] = Xt.t
    df = df.set_index("transfer_n")
    df = df.stack().reset_index()
    df.columns = ["transfer_n", "phenotype", "fraction"]
    df["strategy"] = strategy
    df["method"] = method
    return df

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx