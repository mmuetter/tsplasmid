import pandas as pd
import os
from model_scripts.instrutcion_set_class import InstructionSet
from model_scripts.simulation_class import Simulation
import seaborn as sns


class Simulation_bundle:
    def __init__(self, date, clean=False, location="general_obj"):
        self.instructions = InstructionSet(date)
        self.location = location
        self.clean = clean
        self.pathes = self.instructions.pathes
        self.date = date

    def validation(self, n_repitions):
        simulation = Simulation(
            self.instructions, location=self.location, clean=self.clean)
        dfs = []
        for rep in range(n_repitions):
            print("rep: ", rep)
            simulation.run_simulation()
            df = simulation.timeplot_data
            df["rep"] = rep
            dfs.append(df)
        self.timeplot_results = pd.concat(dfs)

    def variation(self, n_variations):
        dfs = []
        for n_var in range(n_variations):
            self.instructions.randomize_instructions()
            simulation = Simulation(
                self.instructions, location=self.location, clean=self.clean)
            print("var: ", n_var)
            simulation.run_simulation()
            df = simulation.timeplot_data
            df["rep"] = n_var
            dfs.append(df)
        self.timeplot_results = pd.concat(dfs)

    def quickplot(self):
        sns.relplot(
            data=self.timeplot_results,
            x="transfer_n",
            y="f",
            col="strategy",
            col_wrap=3,
            hue="x_hat",
            kind="line"
        )

    def save_timeplot(self, name):
        path = os.path.join(self.pathes["obj"], "simulations")
        if os.path.exists(path) == False:
            os.mkdir(path)
        self.timeplot_results.to_pickle(
            os.path.join(path,  name))
