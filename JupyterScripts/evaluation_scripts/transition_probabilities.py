import numpy as np
import pandas as pd
from itertools import product
import os
import random


def get_transition_probabilities(data, strains, pathes, case: str = "_t_", path_spec: str = "obj", clean=False):
    results = summarize_data(data, strains)
    res = make_weighted_averages(results, strains)
    path = os.path.join(pathes[path_spec], "M"+case)
    write_M(res, strains, path, clean)

    res.to_csv(os.path.join(
        pathes[path_spec], "tranistion_probs"+case+".csv"))
    return res


def write_M(res, X, path, clean):
    for antibiotic in res.antibiotic.unique():
        M = pd.DataFrame()

        for x_hat, x in product(*[X, X]):
            m = res.loc[(res.x_hat == x_hat) & (
                res.x == x) & (res.antibiotic == antibiotic), "weighted_mean"].values[0]
            if type(m) == str:
                M.loc[x_hat, x] = "?"
            else:
                M.loc[x_hat, x] = float(m)

        if clean:
            fpath = os.path.join(path+antibiotic+"_clean.pkl")
            M_clean = clean_M(M, X)
            M_clean.to_pickle(fpath)
        else:
            fpath = os.path.join(path+antibiotic+".pkl")
            M.to_pickle(fpath)


def clean_M(M_raw, strains):
    forbidden_transistions = {
        "U": ["S", "A_r", "B_r", "A&B", "AB_r"],
        "S": ["A_r", "B_r", "A&B", "AB_r"],
        "A_r": ["B_r", "A&B", "AB_r"],
        "B_r": ["A_r", "A&B", "AB_r"],
        "A&B": [],
        "AB_r": []}
    M = M_raw.copy()
    for x in strains:
        if M[x].isin(["?"]).any() == False:
            forbidden = forbidden_transistions[x]
            allowed = list(set(strains) - set(forbidden))
            M.loc[allowed, x] = M.loc[allowed, x]/M.loc[allowed, x].sum()
            M.loc[forbidden, x] = 0
    return M


def summarize_data(data: pd.DataFrame, strains: list):
    data = data[data.x.isin(strains) & data.phenotype.isin(strains)]
    antibiotics = data.treatment_with.unique()
    results = pd.DataFrame()
    t = data.transfer_n.unique()
    plates = data.plate.unique()
    experiments = data.exp.unique()

    for x_in, antibiotic, t, plate, exp in product(*[strains, antibiotics, t, plates, experiments]):
        sub = data[(data.x == x_in) & (data.treatment_with == antibiotic) & (
            data.transfer_n == t) & (data.plate == plate) & (data.exp == exp)]
        if len(sub) > 0:
            for x_hat in strains:
                samples = np.array(sub.phenotype == x_hat)
                row = {"transfer_n": t, "plate": plate, "antibiotic": antibiotic, "exp": exp, "strain_x": x_in,
                       "strain_x_hat": x_hat, "mean": samples.mean(), "weight": len(sub)}
                results = results.append(row, ignore_index=True)
    return results


def make_weighted_averages(results, strains):
    results["mw"] = results["mean"] * results["weight"]
    rows = []
    antibiotics = results.antibiotic.unique()
    for x, x_hat, antibiotic in product(*[strains, strains, antibiotics]):
        sub = results[(results.strain_x == x) & (
            results.strain_x_hat == x_hat) & (results.antibiotic == antibiotic)].copy()
        if len(sub) > 0:
            weight = sub.weight.sum()
            avr = sub.mw.sum()/weight
        else:
            avr = "?"
        row = {"antibiotic": antibiotic, "x": x, "x_hat": x_hat,
               "weighted_mean": avr}
        rows.append(row)
    res = pd.DataFrame.from_records(rows)
    res = res.sort_values(["antibiotic", "x"])
    return res


def get_added_phenotype(row):
    # receives precleaned rows of a dataframe.
    # precleaned: only wells that received a single input drop (no infection)
    if row["transferred_phenotype"]:
        added_strain = row["transferred_phenotype"]
    elif row["turnover_strain"]:
        added_strain = row["turnover_strain_real"]
    return added_strain


def get_clearance_frame(df):
    interesting_columns = ["added_phenotype", "treatment_with", "strategy"]
    df_clearance = df[interesting_columns].copy()
    df_clearance["cleared"] = df.phenotype == "U"
    combis_for_clearance = {
        "S": ["A", "B", "AB"],
        "A_r": ["B", "AB"],
        "B_r": ["A", "AB"]
    }
    df_clearance["n"] = False
    for s in combis_for_clearance:
        for antibiotic in combis_for_clearance[s]:
            df_clearance["n"] = df_clearance["n"] + \
                ((df_clearance.treatment_with == antibiotic)
                 & (df.added_phenotype == s))
    df_clearance = df_clearance[df_clearance["n"]]
    return df_clearance


def get_clearance_props(df_clearance):
    df_total = df_clearance[["added_phenotype", "treatment_with", "n"]].groupby(
        ["added_phenotype", "treatment_with"]).sum().reset_index()
    clearance_props = df_clearance[["added_phenotype", "treatment_with", "cleared"]].groupby(
        ["added_phenotype", "treatment_with"]).sum().reset_index()
    clearance_props["n"] = df_total["n"]
    clearance_props["c_av"] = clearance_props["cleared"]/clearance_props["n"]
    return clearance_props


def round_transition_matrix(df, n=2):
    df_rounded = df
    for col in df.columns:
        adjust = 1
        one = 0
        if df[col].dtype == "float64":
            while not one:
                c_rounded = (adjust * df[col]).round(n)
                one = c_rounded.sum() == 1
                if c_rounded.sum() < 1:
                    adjust += random.random()/100
                elif c_rounded.sum() > 1:
                    adjust -= random.random()/100
            df_rounded[col] = c_rounded
        else:
            df_rounded[col] = df[col]
    return df_rounded
