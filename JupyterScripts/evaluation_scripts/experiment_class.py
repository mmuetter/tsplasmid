import pandas as pd
import numpy as np
from evaluation_scripts.base import build_dictionaries_3, get_pathes_2, rename_strings, df_to_dict
import os
import json
import matplotlib.pyplot as plt


class Experiment:
    def __init__(self, experiment, src="01b_Data", exclude=True):
        self.__dict__.update(build_dictionaries_3())
        self.strategies = list(self.plates.values())
        self.antibiotics = ["none", "A", "B", "AB"]

        if type(experiment) == str:
            self.single = True
            self.pathes = get_pathes_2(experiment+"_output")
            self.exp = experiment
            df = pd.read_pickle(self.pathes["obj"]+os.sep + src + ".pkl")
            df = rename_strings(df, self.rename)
            df["exp"] = experiment
            self.T = df.transfer_n.unique()
            self.T.sort()
        else:
            dfs = []
            for exp in experiment:
                self.single = False
                self.pathes = get_pathes_2(exp+"_output")
                df = pd.read_pickle(self.pathes["obj"]+os.sep + src + ".pkl")
                df = rename_strings(df, self.rename)
                df["exp"] = exp
                dfs.append(df)
            df = pd.concat(dfs)
        self.data_raw = df.copy()
        if exclude:
            self.data = df[df.exclude != True].copy()
        else:
            self.data = df.copy()

        with open(os.path.join(self.pathes["exp"], "exp_pars.json")) as json_data:
            self.exp_pars = json.load(json_data)
        self.strains = ["U", "S", "A_r", "B_r", "A&B", "AB_r"]
        self.get_community()
        self.get_real_x0()
        self.get_community_distribution()
        self.control_wells = self.exp_pars["format"][0]["rwells_bl"]

    def get_community(self):
        community = {}
        total = 0
        for c in self.exp_pars["community"]:
            strain = c["strain"]
            if strain == "wt":
                strain = "S"
            elif strain == "UI":
                strain = "U"
            n = c["n"]
            total += n
            community.update({strain: n})
        for key in community:
            community.update({key: community[key]})
        community.update({"A&B": 0})
        self.community_numbers = community

    def get_community_distribution(self):
        total = np.array(list(self.community_numbers.values())).sum()
        community_distribution = []
        for strain in self.strains:
            community_distribution.append(self.community_numbers[strain]/total)
        self.community_distribution = dict(
            zip(self.strains, community_distribution))

    def get_real_x0(self):
        Data = self.data
        df = Data[Data.turnover_strain.isnull() == False].copy()
        df = df[(df.transfer_n == 0) & (df.plate == 1) & (df.rep == 1)]
        df["x0"] = True
        cols = ["turnover_strain", "x0"]
        x0 = df[cols].groupby(cols[:-1]).count()
        self.x0 = df_to_dict(x0)

    def save_figure(self, name):
        name = self.exp + "_" + name
        file_path = os.path.join(self.pathes["analysis"], name)
        plt.savefig(file_path, bbox_inches='tight', dpi=300)

    def pickle_df(self, df, name):
        file_path = os.path.join(self.pathes["obj"], name)
        df.to_pickle(file_path)
