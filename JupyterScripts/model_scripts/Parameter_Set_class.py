import random
from evaluation_scripts.experiment_class import Experiment
import numpy as np


class Parameter_set:
    def __init__(self, date):
        exp = Experiment(date, exclude=False)
        self.date = date
        self.strains = exp.strains
        self.pathes = exp.pathes
        self.exp_pars = exp.exp_pars
        self.turnover = exp.exp_pars["rates"][0]["turnover"]
        self.infection = exp.exp_pars["rates"][0]["infection"]
        self.community_distribution = exp.community_distribution
        self.exp = exp
        self.control_wells = self.exp_pars["format"][0]["rwells_bl"]

    def randomize_parameter(self, preexisting:bool = True):
        self.turnover  = random.choice([i/100 for i in range(5, 95)])
        self.infection = random.choice([i/100 for i in range(5, 95)])

        community = [random.random() for _ in range(6)]
        community[4] = 0
        if not preexisting:
            community[5] = 0
        community = np.floor(community/np.sum(community)*100)/100
        community[0] = 1 - np.sum(community[1:])
        self.community_distribution = dict(
            zip(self.strains, np.round(community, 2)))

    def __repr__(self):
        return f"Parameter_set(turnover={self.turnover}, infection={self.infection}, community={self.community_distribution})"

