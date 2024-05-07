import seaborn as sns
import pandas as pd
import numpy as np
import os


class Simulation:
    def __init__(self, instruction_set, clean=False, location="general_obj"):
        self.strains = instruction_set.strains
        self.instructions = instruction_set.instructions
        self.simulation_data = self.instructions.copy()
        self.pathes = instruction_set.pathes
        self.M_type = location
        if location != "general_obj":
            name_spec = "_hyb"
        else:
            name_spec = ""
        self.M = self.get_Ms("M"+name_spec+"_t", clean=clean)
        self.M1 = self.get_Ms("M"+name_spec+"_t1", clean=clean)
        self.T = instruction_set.T

    def run_simulation(self):
        #  ADD THE TURNOVER TO ADDED PHENO LIST
        transfer_mask = self.simulation_data.turnover_strain.isnull()
        turnover_mask = transfer_mask == False
        df = self.simulation_data[turnover_mask].copy()
        df["added_phenos"] = df.apply(lambda x: [x.turnover_strain], axis=1)
        self.simulation_data.loc[turnover_mask,
                                 "added_phenos"] = df["added_phenos"]
        # initialize the added pheno list without turnover
        self.simulation_data.loc[transfer_mask, "added_phenos"] = self.simulation_data[transfer_mask].apply(
            lambda x: [], axis=1)

        for t in self.T:
            idx = self.simulation_data.transfer_n == t
            self.df_t = self.simulation_data[idx].copy()

            # ADD TRANSFER Phenos
            self.add_transfer_phenos()

            # ADD INFECTION PHENOS
            self.df_t["added_phenos"] = self.df_t.apply(
                lambda row: self.get_infection_phenos(row), axis=1)

            if t == 1:
                M = self.M1
            else:
                M = self.M

            # get phenotype before treatment
            self.df_t["x"] = self.df_t.apply(
                lambda row: self.get_x(row), axis=1)

            #  get phenotype after treatment
            self.df_t["x_hat"] = self.df_t.apply(
                lambda x: calc_x_hat(x, M, self.strains), axis=1)
            self.simulation_data = merge_sub_df(
                self.simulation_data, self.df_t)
        self.make_timplot_summary()

    def get_Ms(self, name, clean=False):
        path = self.pathes[self.M_type]
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

    def calc_x_mix(self, a, b):
        xa = self.strains.index(a)
        xb = self.strains.index(b)
        infection_matrix = np.array([
            ["U",   "S",    "A_r",  "B_r",  "A&B",  "AB_r"],
            ["S",    "S",    "A_r",  "B_r",  "A&B",  "AB_r"],
            ["A_r",  "A_r",  "A_r",  "A&B",  "A&B",  "AB_r"],
            ["B_r",  "B_r",  "A&B",  "B_r",  "A&B",  "AB_r"],
            ["A&B",  "A&B",  "A&B",  "A&B",  "A&B",  "AB_r"],
            ["AB_r", "AB_r", "AB_r", "AB_r", "AB_r", "AB_r"],
        ])
        return infection_matrix[xa][xb]

    def get_x(self, row):
        x_new = "U"
        for pheno in row["added_phenos"]:
            x_new = self.calc_x_mix(x_new, pheno)
        return x_new

    def get_transfer_phenos(self, row):
        transfer_id = "t"+str(int(row.transfer_n)-1) + "_p" + \
            str(int(row.plate)) + "_rw" + str(int(row.rwell))
        transfer_pheno = self.simulation_data.loc[transfer_id, "x_hat"]
        row.added_phenos.append(transfer_pheno)

    def add_transfer_phenos(self):
        df_t = self.df_t
        transfer_mask = (df_t.transferred_from_well.isnull() == False)
        df = df_t[transfer_mask]
        if not df.empty:
            df.apply(lambda row: self.get_transfer_phenos(row), axis=1)
            self.df_t[transfer_mask] = df

    def get_infection_phenos(self, row):
        added_phenos = row.added_phenos
        for well in row.infected_by:
            added_phenos.append(self.simulation_data.loc[well, "x_hat"])
        return added_phenos

    def quickplot(self):
        sns.relplot(
            data=self.timeplot_data,
            x="transfer_n",
            y="f",
            col="strategy",
            col_wrap=3,
            hue="x_hat",
            kind="line")

    def make_timplot_summary(self):
        df = self.simulation_data
        self.patient_num = patient_num = len(df.rwell.unique())
        df["n"] = True
        cols = ["transfer_n", "strategy", "x_hat", "n"]
        summary = df[cols].groupby(
            cols[:-1]).count().unstack(fill_value=0).stack()
        summary["f"] = summary["n"]/patient_num
        self.timeplot_data = summary.reset_index()


def calc_x_hat(row, M, strains):
    dist = M[row.treatment_with][:, strains.index(row.x)]
    x_hat = np.random.choice(strains, p=list(dist))
    return x_hat


def merge_sub_df(df, df_sub):
    extra_columns = list(set(df_sub.columns) - set(df.columns))
    df[extra_columns] = None
    df.update(df_sub)
    return df
