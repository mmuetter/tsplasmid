import pandas as pd
import numpy as np
import json
import re
import os
from evaluation_scripts.base import FileNames
from evaluation_scripts.well_ident import build_dictionaries


def phenotyping(Data, path):
    _, strategies, _, agar_strain_dict = build_dictionaries()

    agar_file_names = FileNames(".json", path, "pickolo"+os.sep+"json")

    #  Read agar plates
    Data = agar_to_pheno(Data, agar_file_names)

    #  Calculate the phenotype based on the agar-strigns
    agar_strain_dict.update({"?": "?"})
    Data['pheno_num'] = Data.apply(
        lambda x:  strain_arr_to_num(x["pheno_str"]), axis=1)
    Data["phenotype"] = Data.apply(lambda x:  strain_num_to_pheno(
        agar_strain_dict, x["pheno_num"]), axis=1)

    # Label strategy according to plate
    Data['strategy'] = Data.apply(
        lambda x:  strategies["P"+str(x["plate"])], axis=1)

    #  Save strain_dict
    #strain_dict = open(analysis_path+"strain_dict.pkl", "wb")
    #pickle.dump(agar_strain_dict, strain_dict)
    # strain_dict.close()
    return Data


def DisentangleName(name):
    #  384 Well Plates
    if "AI" in name or "BI" in name:
        l1 = name.split("_")
        t = re.findall("\d+", l1[0])
        l2 = l1[-1].split(".")
        barcode = l2[0]
        plate_num = l1[1]

        if "AI" in l1:
            m = "AI"
        elif "BI" in l1:
            m = "BI"
        else:
            print("Something is wrong")
        return int(t[0]), barcode, plate_num, m
    #  Agar Plates
    else:
        l1 = name.split("_")
        t = re.findall("\d+", l1[0])
        l2 = l1[-1].split(".")
        id_num = l2[0]
        plate_num = l1[1]
        treat = l1[2]
        return t[0], plate_num, id_num, treat


def get_transfer_strain(row, Data):
    # Nan is a float -> No infection/transfer
    to_well_id = row["transfer_to_well_id"]
    if type(to_well_id) == str:
        if row.exclude:
            Data.loc[to_well_id, "received_transfer_strain"] = "?"
        else:
            Data.loc[to_well_id,
                     "received_transfer_strain"] = Data.loc[row.name, "phenotype"]
        Data.loc[to_well_id, "contaminated"] = Data.loc[to_well_id,
                                                        "contaminated"] | row.contaminated


def get_infection_strain(row, Data):
    to_well_id = row["infection_to_well_id"]
    if row.exclude:
        strain = "?"
    else:
        strain = row.phenotype
    # Nan is a float -> No infection/transfer
    if (type(to_well_id) == str) & (to_well_id in Data.index):
        # Check if the well was already infected by another strain
        # No-> start list
        if Data.loc[to_well_id, "infected_by_wells"] == None:
            Data.at[to_well_id, "infected_by_wells"] = [row.name]
            Data.at[to_well_id, "infected_by_strains"] = [strain]
        # Else: Append to list
        else:
            Data.at[to_well_id, "infected_by_wells"].append(row.name)
            Data.at[to_well_id, "infected_by_strains"].append(strain)
        Data.loc[to_well_id, "contaminated"] = Data.loc[to_well_id,
                                                        "contaminated"] | row.contaminated


def strain_arr_to_num(strain_arr):
    if (strain_arr is None):
        strain_num = "?"
    elif len(strain_arr) < 4:
        strain_num = "?"
    else:
        strain_arr = "".join(np.array(strain_arr).astype(str))
        if strain_arr == "????":
            strain_num = "?"
        else:
            strain_num = int(strain_arr, 2)
    return strain_num


def strain_num_to_pheno(agar_strain_dict, key):
    if key in agar_strain_dict.keys():
        strain = agar_strain_dict[key]
    else:
        strain = "Fishy"
    return strain


def get_turnover_strain(row, Data, strainplate):
    turnover_id = row["turnover_id"]
    # Nan is a float -> No turnover
    if type(turnover_id) == str:
        Data.loc[row.name,
                 "turnover_strain_real"] = strainplate.loc[turnover_id, "strain_real"]
        Data.loc[row.name, "contaminated"] = strainplate.loc[turnover_id,
                                                             "contaminated"] | Data.loc[row.name, "contaminated"]


def agar_str(ab, a, b, n):
    if (type(ab) == bool) & (type(a) == bool) & (type(b) == bool) & (type(n) == bool):
        return [str(int(ab)), str(int(a)), str(int(b)), str(int(n))]
    else:
        return ["?", "?", "?", "?"]


def agar_to_pheno(Data, agar_file_names):
    # Load plate files
    for index, file in agar_file_names.iterrows():
        t, p, id_num, ab = DisentangleName(file["name"])

        # Load single json file == one agar plate picture (e.g. t = 0, P3, N):
        with open(file["path"]) as f:
            agar_plate = json.load(f)
            agar_plate = agar_plate["growth"]

        # Convert to dataframe with standard IDs
        agar_plate_frame = pd.DataFrame(agar_plate)
        w = "_" + \
            agar_plate_frame["row"].astype(
                str) + agar_plate_frame["col"].astype(str)
        id_col = "t"+t+"_"+p+w
        agar_plate_frame["id"] = id_col
        agar_plate_frame = agar_plate_frame.set_index("id")

        #  If algorithm automatically chose well use that. If corrected manually use that.
        manual = agar_plate_frame.manual != "NA"
        automatic = agar_plate_frame.manual == "NA"
        agar_plate_frame.loc[manual,
                             "real"] = agar_plate_frame.loc[manual.index, "manual"]
        agar_plate_frame.loc[automatic,
                             "real"] = agar_plate_frame.loc[automatic.index, "auto"]

        # Paste Infos to columns in the Data-Dataframe
        Data.loc[agar_plate_frame.index, "agar_"+ab] = agar_plate_frame.real

        # Exclude wells that don't contain any data
        exclude = agar_plate_frame[agar_plate_frame.real == "NA"]
        Data.loc[exclude.index, "exclude"] = True
        Data.loc[exclude.index, "comment"] = "plate unreadable"

    # Convert plate columns to phenostring
    Data['pheno_str'] = Data.apply(lambda x:  agar_str(
        x["agar_AB"], x["agar_A"], x["agar_B"], x["agar_N"]), axis=1).apply(np.array)
    Data = Data.drop(columns=['agar_N', 'agar_B', 'agar_A', 'agar_AB'])
    return Data


def get_x(x_hat):
    if "AB_r" in x_hat:
        x = "AB_r"
    elif "Fishy" in x_hat:
        x = "Fishy"
    elif ("A&B" in x_hat) | (("A_r" in x_hat) & ("B_r" in x_hat)):
        x = "A&B"
    elif "A_r" in x_hat:
        x = "A_r"
    elif "B_r" in x_hat:
        x = "B_r"
    elif "S" in x_hat:
        x = "S"
    elif "?" in x_hat:
        x = "?"
    else:
        x = "UI"
    return x
