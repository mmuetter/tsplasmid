import re
import pandas as pd


def makeID(t, p, row, col):
    id_ = "t"+str(t)+"_P"+str(p)+"_"+str(row)+str(col)
    return id_


def Label_Replicates(Data, plate_rows):
    row_indices = pd.Series(plate_rows).reset_index().set_index(0)['index']
    rep = (row_indices.loc[Data['row']].values %
           2) + ((Data['col'] - 1) % 2) * 2
    Data['rep'] = rep
    return Data


def back_date_id(id_):
    splitted = id_.split("_")
    t_str = splitted[0]

    t = int(re.sub('\D', '', t_str))

    splitted[0] = "t" + str(t-1)
    id_old = "_".join(splitted)
    return id_old


def build_dictionaries():
    # Build dictionary to assign colors
    hue_infos = {"U": "gray", "S": "green", "A": "blue", "B": "orange", "AB": "red", "A&B": "pink", "Fishy": "black",
                 "Fishy_OD": "Brown", "B_r": "orange", "AB_r": "red", "UI": "gray", "?": "white", "WT": "green", "A_r": "blue"}
    #  Build dictionary to assign strategies
    strategies = {"P1": "No treatment", "P2": "Mono A",
                  "P3": "Mono B", "P4": "Combo", "P5": "Cycling", "P6": "Mixing"}
    # Agar Strings to Strains... (AB, A, B, N) -> Binary String -> Decimal -> Strain
    agar_strain_dict = {0: "UI", 5: "A_r",
                        3: "B_r", 15: "AB_r", 7: "A&B", 1: "WT"}
    #  Well Dictionary
    well = 1
    well_dict = {}
    plate_rows = ["A", "B", "C", "D", "E", "F", "G",
                  "H", "I", "J", "K", "L", "M", "N", "O", "P"]
    for row in plate_rows:
        for col in range(1, 25):
            well_dict.update({well: str(row)+str(col)})
            well += 1
    return hue_infos, strategies, plate_rows, agar_strain_dict
