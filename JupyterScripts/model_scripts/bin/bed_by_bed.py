from evaluation_scripts.base import get_pathes_2, build_dictionaries_2
import pandas as pd
import numpy as np
import os
import json
import random


def make_sim(Data, strainplate, strains, sim_number, pathes, location, spec="", clean=False):
    results = pd.DataFrame()
    for i in range(sim_number):
        print("i", i)
        for t in Data.transfer_n.unique():
            if t == 1:
                M = get_Ms("M"+spec+"_t1", pathes[location], clean=clean)
            else:
                M = get_Ms("M"+spec+"_t", pathes[location], clean=clean)
            df = calc_timepoint_val(
                Data, t, M, strainplate, strains).copy()
        df["r"] = df.rep.astype(str) + "_" + str(i)
        df["n"] = True
        results = results.append(df, ignore_index=True)
    return results


def make_sim_altern(Data, strains, sim_number, pathes, location, exp_pars, community, spec="", clean=False, start_at=0):
    results = pd.DataFrame()
    for i in range(start_at, sim_number):
        print("\ni", i)
        print("randomize decisions")
        Data = randomize_decisions(Data, community, exp_pars)
        print("calculate biological stoch")
        for t in Data.transfer_n.unique():
            if t == 1:
                M = get_Ms("M"+spec+"_t1", pathes[location], clean=clean)
            else:
                M = get_Ms("M"+spec+"_t", pathes[location], clean=clean)
            df = calc_timepoint_altern(Data, t, M, strains).copy()
        df["r"] = i
        df["n"] = True
        results = results.append(df, ignore_index=True)
    return results


def randomize_decisions(df, community, exp_pars):
    df["turnover"] = df.apply(lambda x: assign_turnovers(
        x, community, exp_pars["rates"][0]["turnover"]), axis=1)
    df["infected_by"] = df.apply(lambda x: [], axis=1)
    df.apply(lambda x: assign_infected_by(
        x, exp_pars["rates"][0]["infection"], df), axis=1)
    df["plate_well_id"] = "rw" + df.rwell.astype(str) + "_" + df.strategy
    idx = df[df.turnover.isnull()].index
    df.loc[idx, "transfer_well"] = df.loc[idx, :].apply(
        lambda x: get_transfer_well(x, df), axis=1)
    return df


def load_info(exp):
    pathes = get_pathes_2(exp+"_output")
    dictionaries = build_dictionaries_2()
    Data = pd.read_pickle(pathes["obj"]+os.sep+"01b_Data.pkl")
    t_end = max(Data["transfer_n"])+1
    with open(os.path.join(pathes["exp"], "exp_pars.json")) as json_data:
        exp_pars = json.load(json_data)
    strains = ["UI", "S", "A_r", "B_r", "A&B", "AB_r"]
    community, x0 = get_community(exp_pars, strains, pathes)
    return community, x0, strains, exp_pars, dictionaries, pathes, Data, t_end


def calc_x_hat(x, M, strains):
    dist = M[x.treatment_with][:, strains.index(x.x_sim)]
    x_hat = np.random.choice(strains, p=list(dist))
    return x_hat


def calc_timepoint_val(df, t, M, strainplate, strains):
    idx = (df.transfer_n == t)
    df.loc[idx, "x_sim"] = df.loc[idx, ].apply(
        lambda row: get_x_new(row, strains, df, strainplate), axis=1)
    df.loc[idx, "x_hat_sim"] = df.loc[idx, :].apply(
        lambda x: calc_x_hat(x, M, strains), axis=1)
    return df


def get_x_hat_well(well, strainplate, Data):
    if well.split("_")[1] == "S":
        x_hat = strainplate.loc[well, "strain_real"]
        if x_hat == "wt":
            x_hat = "S"
        elif x_hat == "AB":
            x_hat = "AB_r"
    else:
        x_hat = Data.loc[well, "x_hat_sim"]
    return x_hat


def calc_x_mix(a, b, strains):
    xa = strains.index(a)
    xb = strains.index(b)
    infection_matrix = np.array([
        ["UI",   "S",    "A_r",  "B_r",  "A&B",  "AB_r"],
        ["S",    "S",    "A_r",  "B_r",  "A&B",  "AB_r"],
        ["A_r",  "A_r",  "A_r",  "A&B",  "A&B",  "AB_r"],
        ["B_r",  "B_r",  "A&B",  "B_r",  "A&B",  "AB_r"],
        ["A&B",  "A&B",  "A&B",  "A&B",  "A&B",  "AB_r"],
        ["AB_r", "AB_r", "AB_r", "AB_r", "AB_r", "AB_r"],
    ])
    return infection_matrix[xa][xb]


def get_x_new(row, strains, Data, strainplate):
    x_new = "UI"
    for well in row["added_wells"]:
        x_hat = get_x_hat_well(well, strainplate, Data)
        x_new = calc_x_mix(x_new, x_hat, strains)
    return x_new


def calc_timepoint_altern(df, t, M, strains):
    idx = (df.transfer_n == t)
    df.loc[idx, "x_sim"] = df.loc[idx, ].apply(
        lambda row: get_x_new_altern(row, strains, df), axis=1)
    df.loc[idx, "x_hat_sim"] = df.loc[idx, :].apply(
        lambda x: calc_x_hat(x, M, strains), axis=1)
    return df


def get_x_new_altern(row, strains, Data):
    if row.turnover in strains:
        x_new = row.turnover
    else:
        x_new = Data.loc[row.transfer_well, "x_hat_sim"]
    for well in row["infected_by"]:
        x_hat = Data.loc[well, "x_hat_sim"]
        x_new = calc_x_mix(x_new, x_hat, strains)
    return x_new


def get_transfer_well(row, df):
    t = row.transfer_n - 1
    w = row.plate_well_id
    transfer_well = df[(df.transfer_n == t) & (df.plate_well_id == w)]
    return transfer_well.index[0]


def assign_turnovers(x, community, tau):
    if x.transfer_n == 0:
        tau = 1

    r = random.random()
    if r < tau:
        turnover = np.random.choice(
            list(community.keys()), p=list(community.values()))
    else:
        turnover = None
    return turnover


def assign_infected_by(x, beta, df):
    t = x.transfer_n
    if t < max(df.transfer_n):
        df_next = df[(df.transfer_n == t+1) & (df.strategy == x.strategy)]

        if (random.random() < beta):
            infect = np.random.choice(df_next.index)
            infected_by = df.loc[infect, "infected_by"]
            infected_by.append(x.name)
            df.at[infect, "infected_by"] = infected_by


def randomize_decisions(data, community, exp_pars):
    df = data
    interesting_columns = ["rwell", "transfer_n", "treatment_with", "strategy"]
    df = data[(data.turnover_strain != 'bl') & (
        data.rep == 1)][interesting_columns].copy()

    df["turnover"] = df.apply(lambda x: assign_turnovers(
        x, community, exp_pars["rates"][0]["turnover"]), axis=1)
    df["infected_by"] = df.apply(lambda x: [], axis=1)
    df.apply(lambda x: assign_infected_by(
        x, exp_pars["rates"][0]["infection"], df), axis=1)
    df["plate_well_id"] = "rw" + df.rwell.astype(str) + "_" + df.strategy
    idx = df[df.turnover.isnull()].index
    df.loc[idx, "transfer_well"] = df.loc[idx, :].apply(
        lambda x: get_transfer_well(x, df), axis=1)
    return df


def get_Ms(name, path, clean=False):
    if clean:
        c = "_clean"
    else:
        c = ""
    M = {}
    for antibiotic in ["none", "A", "B", "AB"]:
        M.update({
            antibiotic: pd.read_pickle(os.path.join(
                path, name+"_"+antibiotic+c+".pkl")).to_numpy()
        })
    return M


def get_community(exp_pars, strains, pathes, val=False):
    community = {}
    total = 0
    for c in exp_pars["community"]:
        strain = c["strain"]
        if strain == "wt":
            strain = "S"
        n = c["n"]
        total += n
        community.update({strain: n})
    for key in community:
        community.update({key: community[key]})
    community.update({"A&B": 0})
    community, x = sort_phenos(community, strains)

    if val:
        x = get_real_x0(pathes, strains)
    return community, x


def sort_phenos(dictionary, order):
    x = []
    c = {}
    for i, o in enumerate(order):
        if o in dictionary.keys():
            x.append(dictionary[o])
            c.update({o: dictionary[o]/94})
        else:
            x.append(0)
            c.update({o: 0})
    return c, x


def get_real_x0(pathes, strains):
    x0_raw = pd.read_pickle(os.path.join(pathes["obj"], "04c_x0.pkl"))
    x0 = []
    for key in strains:
        if key in x0_raw.index:
            x0.append(x0_raw.loc[key, "x0"])
        else:
            x0.append(0)
    return x0
