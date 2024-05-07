import pandas as pd
import numpy as np
import os
import random


class BbB:
    def __init__(self, exp, M_average=True, clean=False, spec=""):
        self.data = exp.data
        self.T = exp.T
        self.strains = exp.strains
        self.pathes = exp.pathes
        self.spec = spec
        self.exp_pars = exp.exp_pars
        self.turnover = exp.exp_pars["rates"][0]["turnover"]
        self.infection = exp.exp_pars["rates"][0]["infection"]
        self.rates = exp.exp_pars["rates"][0]
        if M_average:
            self.location = self.pathes["general_obj"]
        else:
            self.location = self.pathes["obj"]

        self.clean = clean

    def generate_parameter(self, R0, res):
        if R0 < 1:
            tau = random.random()
            beta = R0*tau
        else:
            beta = random.random()
            tau = beta/R0
        self.turnover = tau
        self.infection = beta

        xA = random.random()*res
        xB = res-xA
        xS = (1-res)*random.random()
        xU = 1-res-xS
        X = np.array([xU, xS, xA, xB])
        self.community = dict(zip(self.strains[0:4], X))

    def randomize_decisions(self):
        self.data["turnover"] = self.data.apply(
            lambda x: self.assign_turnovers(x), axis=1)
        self.data["infected_by"] = self.data.apply(lambda x: [], axis=1)
        self.assign_infected_by()
        self.data["plate_well_id"] = "rw" + \
            self.data.rwell.astype(str) + "_" + self.data.strategy
        idx = self.data[self.data.turnover.isnull()].index
        self.data.loc[idx, "transfer_well"] = self.data.loc[idx, :].apply(
            lambda x: get_transfer_well(x, self.data), axis=1)

    def assign_turnovers(self, x):
        if x.transfer_n == 0:
            tau = 1
        else:
            tau = self.rates["turnover"]

        r = random.random()
        if r < tau:
            turnover = np.random.choice(
                list(self.community.keys()), p=list(self.community.values()))
        else:
            turnover = None
        return turnover

    def make_sim_altern(self, sim_number, spec="", clean=False, start_at=0):
        results = []
        for i in range(start_at, sim_number):
            print("calculate biological stoch")
            for t in self.data.transfer_n.unique():
                if t == 1:
                    M = get_Ms("M"+spec+"_t1", self.location, clean=clean)
                else:
                    M = get_Ms("M"+spec+"_t", self.location, clean=clean)
                df = calc_timepoint_altern(
                    self.Data, t, M, self.strains).copy()
            df["r"] = i
            df["n"] = True
            results.append(df)
        results = pd.DataFrame().from_records(results)
        self.sim_results = results

    def assign_infected_by(self):
        self.data.apply(lambda x: self.assign_infected_by_i(x), axis=1)

    def assign_infected_by_i(self, x):
        t = x.transfer_n
        if t < max(self.data.transfer_n):
            df_next = self.data[(self.data.transfer_n == t+1) &
                                (self.data.strategy == x.strategy)]

            if (random.random() < self.rates["infection"]) & (df_next.empty == False):
                infect = np.random.choice(df_next.index)
                infected_by = self.data.loc[infect, "infected_by"]
                infected_by.append(x.name)
                self.data.at[infect, "infected_by"] = infected_by


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
