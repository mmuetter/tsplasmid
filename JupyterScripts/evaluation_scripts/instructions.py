import pandas as pd
import numpy as np
from evaluation_scripts.well_ident import back_date_id
import math


def Import_Raw_Instructions(path, instruction_file_names):
    # read files in the list and concatenate them
    dfs = [pd.read_csv(p) for p in instruction_file_names["path"]]
    instructions = pd.concat(dfs, ignore_index=True)

    # Create the ID column
    t = 't' + instructions["transfer_n"].astype(str)
    p = "_P" + instructions["plate"].astype(str)
    w = "_" + \
        instructions["row"].astype(str).str.cat(
            instructions["col"].astype(str))
    id_ = t.str.cat([p, w], sep='')

    instructions["ID"] = id_
    instructions = instructions.set_index("ID")

    return instructions


def get_turnover_id(strainplate, row):
    # Check if turnover happens
    if type(row["turnover_strain"]) == str:
        well = strainplate[(strainplate["transfer"] == row["transfer_n"]) & (
            strainplate["rep"] == row["rep"]) & (strainplate["strain_wanted"] == row["turnover_strain"])]
        id_ = well.index[0]
    else:
        id_ = None
    return id_


def r_well_to_id_infection(Data):
    t_end = max(Data.transfer_n)
    data = Data[(Data.transfer_n > 0)]
    for i, row in data.iterrows():
        # Infection R_well to ID
        #  Only if this well infects someone
        if np.isnan(row.infection_to_well) != True:
            infected_well = Data[(Data["rwell"] == row["infection_to_well"]) & (
                Data["transfer_n"] == row["transfer_n"]) & (Data["rep"] == row["rep"]) & (Data["plate"] == row["plate"])]
            Data.loc[back_date_id(
                i), "infection_to_well_id"] = infected_well.index
    return Data


def infection_to_well_id(row, Data):
    if math.isnan(row.infection_to_well) == False:
        t_i = row.transfer_n + 1
        r_w = row.infection_to_well
        r = row.rep
        p = row.plate
        df = Data[(Data.rwell == r_w) & (Data.rep == r) & (Data.plate == p)]
        well = df.row.unique()[0] + str(df.col.unique()[0])
        id_ = "t" + str(t_i) + "_P" + str(p) + "_" + well
        return id_
    else:
        return float("nan")


def get_transfer_id(x):
    # t0 -> Initializing -> No transfer
    transfer_id = None
    if x.transfer_n > 0:
        # No Turnover -> Transfer
        if type(x.turnover_strain) != str:
            #print(x.name, "->", ef.back_date_id(x.name))
            transfer_id = back_date_id(x.name)
    return transfer_id


def transfer_from_to_id(row, Data):
    from_id = row["transfer_from_well_id"]
    if type(from_id) == str:
        Data.loc[from_id, "transfer_to_well_id"] = row.name
