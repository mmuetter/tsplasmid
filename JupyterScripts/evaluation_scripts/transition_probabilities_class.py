import numpy as np
import pandas as pd
from itertools import product
import os

from evaluation_scripts.transition_probabilities import make_weighted_averages, clean_M, round_transition_matrix


class Transition_probabilities:
    def __init__(self, experiment, t1=False, write=False):
        self.experiment = experiment
        self.strains = experiment.strains
        self.pathes = experiment.pathes
        self.antibiotics = experiment.antibiotics
        if t1:
            self.case = "_t1"
        else:
            self.case = "_t"

        df = experiment.data[(experiment.data.transfer_n == 1) == t1].copy()
        self.df = df

        self.results = self.summarize_data()
        self.weighted_averages = make_weighted_averages(
            self.results, experiment.strains)

        if self.experiment.single:
            self.path = experiment.pathes["obj"]
        else:
            self.path = experiment.pathes["general_obj"]

        self.get_M(write)

        self.results.to_csv(os.path.join(
            self.path, "tranistion_probs"+self.case+".csv"))

    def summarize_data(self):
        strains = self.experiment.strains
        data = self.df
        data = data[data.x.isin(strains) & data.phenotype.isin(strains)]
        antibiotics = data.treatment_with.unique()
        t = data.transfer_n.unique()
        plates = data.plate.unique()
        experiments = data.exp.unique()
        rows = []
        for x_in, antibiotic, t, plate, exp in product(*[strains, antibiotics, t, plates, experiments]):
            sub = data[(data.x == x_in) & (data.treatment_with == antibiotic) & (
                data.transfer_n == t) & (data.plate == plate) & (data.exp == exp)]
            if len(sub) > 0:
                for x_hat in strains:
                    samples = np.array(sub.phenotype == x_hat)
                    row = {"transfer_n": t, "plate": plate, "antibiotic": antibiotic, "exp": exp, "strain_x": x_in,
                           "strain_x_hat": x_hat, "mean": samples.mean(), "weight": len(sub)}
                    rows.append(row)
        results = pd.DataFrame.from_records(rows)
        return results

    def get_M(self, write=False):
        X = self.strains
        transition_matrices = {}
        for antibiotic, clean in product(self.antibiotics, [True, False]):
            M = pd.DataFrame()
            for x_hat, x in product(*[X, X]):
                m = self.weighted_averages.loc[(self.weighted_averages.x_hat == x_hat) & (
                    self.weighted_averages.x == x) & (self.weighted_averages.antibiotic == antibiotic), "weighted_mean"].values[0]
                if type(m) == str:
                    M.loc[x_hat, x] = "?"
                else:
                    M.loc[x_hat, x] = float(m)

            if clean:
                M_clean = clean_M(M, X)
                transition_matrices.update(
                    {"M" + self.case + "_" + antibiotic+"_clean": round_transition_matrix(M_clean)})
            else:
                transition_matrices.update(
                    {"M" + self.case + "_" + antibiotic: round_transition_matrix(M)})
        self.transition_matrices = transition_matrices
        self.__dict__.update(transition_matrices)
        if write:
            self.write_M()

    def write_M(self):
        for key in self.transition_matrices:
            fpath = os.path.join(self.path, key + ".pkl")
            M = self.transition_matrices[key]
            M.to_pickle(fpath)
