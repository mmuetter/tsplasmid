import pandas as pd
import pickle
import os


class EvalSimGrid:
    def __init__(self, path, parameter_names=["turnover", "infection", "cA", "cB", "cS", "cU", "cAB"], name_add="_ab"):
        self.path = path
        self.file_path = os.path.join(
            path, "sensitivityAnalysis" + name_add + ".pkl")
        self.parameter_names = parameter_names

    def evaluate_parameter_sets(self, simulations, parameter):
        self.simulations = simulations
        self.strategies = simulations[list(simulations.keys())[
            0]].strategy.unique()

        self.parameter = parameter
        self.evaluations = []
        df = []
        for simulation, parameter in zip(list(self.simulations.values()), list(self.parameter.values())):
            parameter_eval = Parameter_Point_Eval(simulation, parameter)
            self.evaluations.append(parameter_eval)
            df.append(parameter_eval.data)
        self.data = pd.concat(df)
        self.strategies = parameter_eval.strategies
        self.data = self.data.rename(columns={
            "U": "cU",
            "S": "cS",
            "A_r": "cA",
            "B_r": "cB",
            "AB_r": "cAB"
        })

    def save(self):
        with open(self.file_path, 'wb') as file:
            pickle.dump(self, file)

    def load(self):
        with open(self.file_path, 'rb') as file:
            obj = pickle.load(file)
        self.__dict__.update(obj.__dict__)


class Parameter_Point_Eval:
    def __init__(self, simulation, parameter_set):
        self.simulation = simulation
        self.parameter_set = parameter_set
        self.summerize_simulation()
        self.add_parameter()

    def summerize_simulation(self, exclude_noTreatment=True):
        self.t_end = t_end = max(self.simulation.transfer_n)
        sim_end = self.simulation[self.simulation.transfer_n.isin(
            [t_end, t_end-1, t_end-2, t_end-3])]
        sim_U = sim_end[sim_end.x_hat == "U"]
        if exclude_noTreatment:
            sim_U = sim_U[sim_U.strategy != "No treatment"]

        self.data = sim_U[["strategy", "f"]].groupby(
            ["strategy"]).mean().reset_index()
        self.strategies = self.data["strategy"].unique()
        self.data = self.data.set_index('strategy').T

    def add_parameter(self):
        for key, value in self.parameter_set.items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    if "winning" not in sub_key:
                        self.data[f'{sub_key}'] = [sub_value]
            else:
                # For regular parameters, add a new column with the parameter name
                self.data[key] = [value]
