import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import re


class well_info:
    def __init__(self, well, track_pheno, df, strainplate, pedigree=False):
        df = df.copy()

        x = df.loc[well, ]

        loss_a_dict = {
            "AB_r": ["B_r", "S"],
            "A_r": ["S"]}
        loss_b_dict = {
            "AB_r": ["A_r", "S"],
            "B_r": ["S"]}

        allowed_mix = ["UI", "S"]
        # If the well is mixed (intended or unintended) with a well that contains other plasmids don't count
        self.count = True
        self.stop = False
        self.comments = []
        self.well_id = "P"+str(x.plate) + "_"+x.row+str(x.col)

        # PEDIGREE
        if bool(pedigree) == False:
            pedigree = [x.turnover_id]

        #  WELL MIXING
        added_wells = x.added_wells.copy()
        added_wells.remove(pedigree[-1])
        if bool(added_wells):
            for well in added_wells:
                if "S" in well:
                    strain = strainplate.loc[well, "strain_real"]
                else:
                    strain = df.loc[well, "phenotype"]

                if ((strain in allowed_mix) == False) & (well not in pedigree):
                    self.count = False
                    self.stop = True
                    self.comments.append(
                        "mixing with other plasmid wells:" + strain)

        #  BASE INFO
        self.well_id = x.name
        self.strategy = x.strategy
        self.track_pheno = track_pheno
        self.x = x
        self.pheno = x.phenotype
        self.transferred_phenotype = x.transferred_phenotype

        # CHECK FOR PLASMID LOSS OR UNINTENDED INFECTION
        self.loss_a = False
        self.loss_b = False
        self.loss_co = False

        if track_pheno != x.phenotype:
            if (track_pheno == "A_r") | (track_pheno == "AB_r"):
                self.loss_a = x.phenotype in loss_a_dict[track_pheno]
                if self.loss_a:
                    self.stop = True
                    self.comments.append("pA lost")
                else:
                    self.stop = True
                    self.count = False
                    self.comments.append("fishy change of pheno")

            if (track_pheno == "B_r") | (track_pheno == "AB_r"):
                self.loss_b = x.phenotype in loss_b_dict[track_pheno]
                if self.loss_b:
                    self.stop = True
                    self.comments.append("pB lost")
                else:
                    self.stop = True
                    self.count = False
                    self.comments.append("fishy change of pheno")

            if (track_pheno == "AB_r"):
                self.loss_co = x.phenotype == "A&B"
                self.comments.append("coexistence lost")
                self.stop = True
                self.count = True
                self.result = "loss co"
                self.comment = "loss of double resistant bacteria"
        # CHILDREN
        self.children = []
        if self.stop == False:
            if bool(x.transfer_to_well_id) & (type(x.transfer_to_well_id) == str):
                self.children.append(x.transfer_to_well_id)
            if bool(x.infection_to_well_id) & (type(x.infection_to_well_id) == str):
                self.children.append(x.infection_to_well_id)

        # UPDATE PEDIGREE
        p = pedigree.copy()
        p.append(x.name)
        self.pedigree = p

        # RESULT
        self.result = "open_end"
        if self.stop & self.count:
            if self.loss_a:
                self.result = "loss_a"
            if self.loss_b:
                self.result = "loss_b"
            if self.loss_co:
                self.result = "loss_co"


def get_children(w, track_pheno, Data, strainplate):
    children = {}
    if bool(w.children):
        for child in w.children:
            children.update(
                {child: well_info(child, track_pheno, Data, strainplate, w.pedigree)})
    return children


def get_new_generation(g, track_pheno, Data, strainplate):
    new_gen = {}
    for g_i in g.values():
        children = g_i.children
        if bool(children):
            for child in children:
                new_gen.update({child: well_info(
                    child, track_pheno, Data, strainplate, pedigree=g_i.pedigree)})
    return new_gen


def get_well_history(start_well, track_pheno, Data, strainplate):
    t = 1
    g = {start_well: well_info(
        start_well, track_pheno, Data.copy(), strainplate)}
    generations = {t: g}

    while bool(g):
        g = get_new_generation(g, track_pheno, Data, strainplate)
        generations.update({t: g})
        t += 1
    return generations


def get_endpoints(generations):
    endpoints = pd.DataFrame(
        columns=["start_well", "transfers", "pedigree", "result", "strategy"])

    for k, g in zip(generations.keys(), generations.values()):
        for l in g.values():
            if l.stop:
                row = {"start_well": l.pedigree[0],
                       "end_well": l.pedigree[-1],
                       "result": l.result,
                       "pedigree": l.pedigree,
                       "transfers": len(l.pedigree)-1,
                       "strategy": l.strategy}
                endpoints = endpoints.append(row, ignore_index=True)
    return endpoints


def show_well_history(endpoints):
    colors = sns.color_palette("pastel")
    reasons = {
        "open_end": "v",
        "loss_a": "x",
        "loss_b": "*",
        "loss_co": "D"
    }

    for i, row in endpoints.iterrows():
        gen = list(range(-1, row.transfers))
        wells = row.pedigree
        plt.plot(gen[:-1], wells[:-1], ".", markersize=15, color=colors[i])
        plt.plot(gen, wells, ":",  color=colors[i])
        plt.plot(gen[-1], wells[-1], reasons[row.result],
                 markersize=15,   color=colors[i])


def show_well_histories(endpoints, figsize=(12, 8)):
    plt.figure(figsize=figsize)

    colors = sns.color_palette("viridis", len(endpoints))
    reasons = {
        "open_end": "v",
        "loss_a": "x",
        "loss_b": "*",
        "loss_co": "D"
    }

    for j, row in endpoints.iterrows():
        i = list(endpoints.index).index(j)
        pedigree = row.pedigree
        T = []
        W = []
        for ped in pedigree[1::]:
            t, p, w = analyse_name(ped)
            T.append(t)
            W.append(w)

        plt.plot(T[:-1], W[:-1], ".", markersize=15, color=colors[i])
        plt.plot(T, W, ":",  color=colors[i])
        plt.plot(T[-1], W[-1], reasons[row.result],
                 markersize=15,   color=colors[i])


def analyse_name(s):
    t, p, w = s.split("_")
    t = int(re.findall(r'\d+', t)[0])
    p = int(re.findall(r'\d+', p)[0])
    return t, p, w
