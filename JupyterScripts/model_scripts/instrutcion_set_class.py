import numpy as np
from model_scripts.Parameter_Set_class import Parameter_set
from evaluation_scripts.base import random_array


class InstructionSet:
    def __init__(self, date=None, parameter=None):
        if parameter:
            self.parameter = parameter
        else:
            self.parameter = parameter = Parameter_set(date)
        self.__dict__.update(parameter.__dict__)
        self.T = self.exp.T
        self.T.sort()
        self.make_infection_distribution()
        # Prep instructions
        self.get_instructions()

    def make_infection_distribution(self):
        infection_distribution = 96*[1/96]
        rwells = list(np.array(range(len(infection_distribution))) + 1)
        for cw in self.control_wells:
            infection_distribution[cw-1] = 0
        infection_distribution = infection_distribution / \
            np.sum(infection_distribution)
        self.infection_distribution = dict(zip(rwells, infection_distribution))

    def randomize_instructions(self):
        #  Randomize infection to and turnover strain
        self.randomize_base_instructions()
        # Handle infections
        self.assign_infected_by()
        #  Get Transferred_from_well
        self.transferred_from_well()

    def randomize_base_instructions(self):
        df = self.instructions.copy()
        #  t0
        mask = df.transfer_n == 0
        length = mask.sum()
        turnovers = random_array(1, self.community_distribution, length)
        infection_to = random_array(
            self.infection, self.infection_distribution, length)
        df.loc[mask, "turnover_strain"] = turnovers
        df.loc[mask, "infection_to_well"] = infection_to
        # t!=0
        mask = df.transfer_n != 0
        length = mask.sum()
        turnovers = random_array(
            self.turnover, self.community_distribution, length)
        infection_to = random_array(
            self.infection, self.infection_distribution, length)
        # save
        df.loc[mask, "turnover_strain"] = turnovers
        df.loc[mask, "infection_to_well"] = infection_to
        self.instructions = df

    def get_instructions(self):
        #  Copy Data
        df = self.parameter.exp.data.copy()
        #  Drop Control Wells
        index = df[df.rwell.isin(self.control_wells)].index
        df = df.drop(index)
        self.instructions = df
        # Set index right
        self.set_index()
        # Reduce to Rwells
        self.condense()
        # Use the real turnover phenos (incl. contamination info)
        self.instructions["turnover_strain"] = self.instructions["turnover_strain_real"]

        # Handle infections
        self.assign_infected_by()
        #  Get Transferred_from_well
        self.transferred_from_well()

        #  Only keep necessary columns
        cols = ["transfer_n", "strategy", "rwell", "infection_to_well", "treatment_with",
                "turnover_strain", "transferred_from_well", "plate", "infected_by"]
        self.instructions = self.instructions[cols]

        # Remove inconsistencies
        self.instructions.turnover_strain = self.instructions.turnover_strain.replace(
            'AB', 'AB_r')
        self.instructions_original = self.instructions.copy()

    def set_index(self):
        df = self.instructions.reset_index()
        df["patient_id"] = "t" + df.transfer_n.astype(
            str) + "_p" + df.plate.astype(str) + "_rw" + df.rwell.astype(str)
        df.set_index("patient_id", drop=True, inplace=True)
        self.instructions = df

    def transferred_from_well(self):
        self.instructions.transferred_from_well = None
        df = self.instructions
        df = df[df.turnover_strain.isnull()]
        self.instructions.loc[list(
            df.index), "transferred_from_well"] = df.rwell

    def condense(self):
        df = self.instructions
        self.instructions = df[df.rep == 0]

    def assign_infected_by(self):
        df = self.instructions
        df["infected_by"] = df.apply(lambda x: [], axis=1)
        mask = df.infection_to_well.isnull() == False
        df[mask].apply(lambda x: self.infected_by(x), axis=1)

    def infected_by(self, x):
        ti = x["transfer_n"]+1
        if ti <= max(self.T):
            infect_idx = "t" + str(ti) + "_p" + \
                str(x["plate"]) + "_rw" + str(int(x["infection_to_well"]))
            self.instructions.loc[infect_idx, "infected_by"].append(x.name)
