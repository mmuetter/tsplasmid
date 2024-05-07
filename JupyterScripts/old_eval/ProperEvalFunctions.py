import pandas as pd
import os 
import re 
import numpy as np
import math
import json
from ast import literal_eval
from random import random

def summerize(df, case):
    summary_1 = df.groupby(["transfer_n", "strategy", "replicate", case]).count()

    strains = summary_1.index.get_level_values(case).unique()
    for t in summary_1.index.get_level_values("transfer_n").unique():
        for s in summary_1.index.get_level_values("strategy").unique():
            for r in summary_1.index.get_level_values("replicate").unique():
                tmp = summary_1.loc[(t, s, r)]
                for strain in strains:
                    if strain not in tmp.index:
                        summary_1.loc[(t, s, r, strain)] = 0

    summary_2 = summary_1.reset_index()
    summary_3 = summary_2.rename(columns={"row": "count", "col" : "fraction"})
    #summary_3 = summary_3[["transfer_n", "strategy", "replicate", "treated_phenotype_machine", "count", "fraction"]]
    for i, line in summary_3.iterrows():
        #summary.loc[i, "transfer_n"] = "t" + str(line["transfer_n"]) 
        summary_3.loc[i, "fraction"] = line["count"]/96*100
    return summary_3

def back_date_id(id_):
    splitted = id_.split("_")
    t_str = splitted[0]

    t = int(re.sub('\D', '', t_str))

    splitted[0] ="t"+ str(t-1)
    id_old = "_".join(splitted)
    return  id_old


def get_transfer_id(x):
    # t0 -> Initializing -> No transfer
    transfer_id = None
    if x.transfer_n > 0:
        # No Turnover -> Transfer
        if type(x.turnover_strain) != str:
            #print(x.name, "->", ef.back_date_id(x.name))
            transfer_id = back_date_id(x.name)
    return transfer_id

def r_well_to_id_infection(Data,strainplate):
    t_end = max(Data.transfer_n)
    data = Data[(Data.transfer_n > 0)]
    for i, row in data.iterrows():
        ## Infection R_well to ID
        ### Only if this well infects someone
        if np.isnan(row.infection_to_well) != True:
            infected_well = Data[(Data["rwell"] == row["infection_to_well"]) & (Data["transfer_n"] == row["transfer_n"]) & (Data["replicate"] == row["replicate"]) & (Data["plate"] == row["plate"])]
            Data.loc[pef.back_date_id(i), "infection_to_well_id"] = infected_well.index
    return Data


def add_pheno_str(a1, a2):
    a_sum = np.array(a1)+np.array(a2)
    a_true = a_sum>0
    return a_true.astype(int)

def superconjugation(strain_arr, probability):

    if list(strain_arr) == [0,1,1,1]:
        if random() < probability:
            strain_arr = [1, 1, 1, 1]
    return strain_arr

from random import *
def treat_phenotype(original_arr, treatment, clearence_rates, strain_dict):
    
    strain = pef.strain_num_to_pheno(strain_dict, pef.strain_arr_to_num(original_arr))   
    if strain == "WT":
        strain = "wt"
       
    ## Calculate the clearence rate for the strain treatment combination
    if strain == "A&B":
        clearence_rate_A = float(clearence_rates[(clearence_rates.treatment_with == treatment) & (clearence_rates.single_input_strain == "A_r")]["clearence_rate"])
        clearence_rate_B = float(clearence_rates[(clearence_rates.treatment_with == treatment) & (clearence_rates.single_input_strain == "B_r")]["clearence_rate"])
        rA = random()
        rB = random()
        cA = clearence_rate_A>rA
        cB = clearence_rate_B>rB
        treatment_arr = [0, cA, cB, cA&cB]
    elif strain == "AB_r":
        try:
            clearence_rate = float(clearence_rates[(clearence_rates.treatment_with == treatment) & (clearence_rates.single_input_strain == strain)]["clearence_rate"])
        except:
            clearence_rate = 0
        rAB = random()
        c_AB = clearence_rate>rAB
        treatment_arr = [c_AB, c_AB, c_AB, c_AB]
    elif (strain == "UI") | (strain =="bl"):
        treatment_arr = [0,0,0,0]
    else:
        #print(treatment, strain)
        clearence_rate = float(clearence_rates[(clearence_rates.treatment_with == treatment) & (clearence_rates.single_input_strain == strain)]["clearence_rate"])
        r = random()
        c = clearence_rate>r
        if strain == "A_r":
            treatment_arr = [0,c,0,c]
        elif strain == "B_r":
            treatment_arr = [0,0,c,c]
        elif strain == "wt":
            treatment_arr = [0,0,0,c]
        else:
            treatment_arr = [0,0,0,0]

    # Create tretead Phenotype     
    ### Check if a "real strain" enters    
    if original_arr is None:
        treated_arr = None
    elif len(original_arr) < 4:
        treated_arr = None
    ## Should be the case
    else:
        treated_arr = np.array(original_arr) - np.array(treatment_arr)
        
        
    return treated_arr


def expect_phenotype_machine(row, Data, strainplate):
    expected_phenotype = row.expected_pheno_str_machine
    
    ## iterate for all input wells
    #print("\n \n Receiving Well:", row.name, "All input wells", row.input_wells, "until now the phenostr is:", expected_phenotype)
    for well in row.input_wells:
        ### check if the id refers to strainplate or assayplate (Data)
        #print("Well", well)
        if "S" in well:
            input_str = strainplate.loc[well, "strain_real_str"]
            #print("Looking for Strainplate:", input_str)
        else:
            input_str = Data.loc[well, "treated_pheno_str_machine"]
            #print("Looking for assay plates", input_str)
            
            ## Check if we have complete information of the input wells... One well with incomplete information leads to unknown results
        #print("Input String:", input_str, "previous_Expected_str", expected_phenotype)
        if len(input_str) + len(expected_phenotype) == 8:
            #print(row.name, expected_phenotype, input_str)
            new_expected_phenotype = pef.add_pheno_str(expected_phenotype, input_str)
            #print("Well",well, ": Expected:", expected_phenotype, "+ Input well", input_str, " = ", new_expected_phenotype)
            expected_phenotype = new_expected_phenotype
        else:
            expected_phenotype = []
        #print(expected_phenotype)
        #print(expected_phenotype)
        Data.at[row.name, "expected_pheno_str_machine"] = expected_phenotype
#Data["expected_pheno_str_machine"] = Data["expected_pheno_str_machine"].apply(lambda x: [0, 0, 0, 0]) 


def expect_phenotype(row, Data, strainplate):
    #print("\n\n", i)
    expected_phenotype = row.expected_phenotype
    
    ## iterate for all input wells
    for well in row.input_wells:
        ### check if the id refers to strainplate or assayplate (Data)
        if "S" in well:
            input_str = strainplate.loc[well, "strain_real_str"]
        else:
            input_str = Data.loc[well, "pheno_str"]
        
        ## Check if we have complete information of the input wells... One well with incomplete information leads to unknown results
        if len(input_str) + len(expected_phenotype) == 8:
            #print(row.name, expected_phenotype, input_str)
            new_expected_phenotype = add_pheno_str(expected_phenotype, input_str)
            #print("Well",well, ": Expected:", expected_phenotype, "+ Input well", input_str, " = ", new_expected_phenotype)
            expected_phenotype = new_expected_phenotype
        else:
            expected_phenotype = []
        Data.at[row.name, "expected_pheno_str"] = expected_phenotype

def create_binary_array(num):
    array = [0, 0, 0, 0]
    if num >= 0:
        b = bin(int(num)).replace("0b", "")
        pos = [3,2,1,0]
        i_pos = 1
        for p in pos:
            if i_pos <= len(b):
                array[p] = int(b[-i_pos])
                i_pos = i_pos + 1
    return array

def strain_to_pheno_num(inv_strain_dict, strain):
    ## If Nan == float do nothing
    if type(strain) == str:
        strain_num =  inv_strain_dict[strain]
    else:
        strain_num = None
    return strain_num

def repair_array(row, Data):
    string = row.pheno_str
    s = "".join(string)
    numbers = []
    for word in s:
        if word.isdigit():
            numbers.append(int(word))
    #print(string, numbers)
    Data.at[row.name, "pheno_str"] = numbers

def input_wells(row, Data):
    input_wells = []
    ## Add turnover wells
    if type(row["turnover_id"]) == str:
        input_wells.append(row["turnover_id"])

    ## Add transfer wells
    if type(row["transfer_from_well_id"]) == str:
        input_wells.append(row["transfer_from_well_id"])
    
    ## Add infection wells
    if type(row["infected_by_wells"]) == str:
        input_wells.extend(literal_eval(row["infected_by_wells"]))
    Data.at[row.name, "input_wells"] = input_wells

def get_turnover_id(strainplate, row):
    ## Check if turnover happens
    if type(row["turnover_strain"]) == str:
        well = strainplate[(strainplate["transfer"] == row["transfer_n"]) & (strainplate["replicate"] == row["replicate"]) & (strainplate["strain_wanted"] == row["turnover_strain"])]
        id_ = well.index[0]
    else:
        id_ = None
    return id_

def get_transfer_strain(row, Data):
    ## Nan is a float -> No infection/transfer
    to_well_id = row["transfer_to_well_id"]
    if type(to_well_id) == str:
        Data.loc[to_well_id, "received_transfer_strain"] = Data.loc[row.name, "phenotype"]
        Data.loc[to_well_id, "contaminated"] = Data.loc[to_well_id, "contaminated"] | Data.loc[row.name, "contaminated"]

def get_infection_strain(row, Data):
    to_well_id = row["infection_to_well_id"]
    strain = row.phenotype
    ## Nan is a float -> No infection/transfer
    if type(to_well_id) == str:
        ## Check if the well was already infected by another strain
        # No-> start list
        if Data.loc[to_well_id, "infected_by_wells"] == None:
            Data.at[to_well_id, "infected_by_wells"] = [row.name]
            Data.at[to_well_id, "infected_by_strains"] = [strain]
        # Else: Append to list
        else:
            Data.at[to_well_id, "infected_by_wells"].append(row.name)
            Data.at[to_well_id, "infected_by_strains"].append(strain)
        Data.loc[to_well_id, "contaminated"] = Data.loc[to_well_id, "contaminated"] | Data.loc[row.name, "contaminated"]

        
def get_turnover_strain(row, Data, strainplate):
    turnover_id = row["turnover_id"]
    ## Nan is a float -> No turnover
    if type(turnover_id) == str:
        Data.loc[row.name, "turnover_strain_real"] = strainplate.loc[turnover_id, "strain_real"]
        Data.loc[row.name, "contaminated"] = strainplate.loc[turnover_id, "contaminated"] | Data.loc[row.name, "contaminated"]

def strain_arr_to_num(strain_arr) :
    if (strain_arr is None):
        strain_num =  "?"
    elif len(strain_arr) < 4:
        strain_num = "?" 
    else:
        strain_arr = "".join(np.array(strain_arr).astype(str))
        if strain_arr == "????":
            strain_num = "?"
        else:
            strain_num = int(strain_arr,2)
    return strain_num

def strain_num_to_pheno(agar_strain_dict, key):
    if key in agar_strain_dict.keys():
        strain = agar_strain_dict[key]
    else:
        strain = "Fishy"
    return strain

def agar_str(ab, a, b, n):
        if (type(ab) == bool) & (type(a) == bool) & (type(b) == bool) & (type(n) == bool):
            return [str(int(ab)), str(int(a)), str(int(b)), str(int(n))]
        else:
            return ["?", "?", "?", "?"]

def Agar_Plate_Files_to_Pheno_strings(Data, agar_file_names):
    plate_string_dict = {"AB" : 0, "A" : 1, "B" : 2, "N" : 3} 
    
    # Load plate files
    for index, file in agar_file_names.iterrows():
        t, p, id_num, ab = pef.DisentangleName(file["name"])

        # Load single json file == one agar plate picture (e.g. t = 0, P3, N): 
        with open(file["path"]) as f:
            agar_plate = json.load(f)
            agar_plate = agar_plate["growth"]

        ## Convert to dataframe with standard IDs
        agar_plate_frame = pd.DataFrame(agar_plate)
        w = "_" + agar_plate_frame["row"].astype(str) + agar_plate_frame["col"].astype(str)
        id_col = "t"+t+"_"+p+w
        agar_plate_frame["id"] = id_col
        agar_plate_frame = agar_plate_frame.set_index("id")

        ## If algorithm automatically chose well use that. If corrected manually use that.
        manual = agar_plate_frame.manual != "NA"
        automatic = agar_plate_frame.manual == "NA"
        agar_plate_frame.loc[manual, "real"] =agar_plate_frame.loc[manual.index,"manual"] 
        agar_plate_frame.loc[automatic, "real"] =agar_plate_frame.loc[automatic.index,"auto"] 

        ## Paste Infos to columns in the Data-Dataframe
        Data.loc[agar_plate_frame.index,"agar_"+ab] = agar_plate_frame.real

        ## Exclude wells that don't contain any data
        exclude = agar_plate_frame[agar_plate_frame.real == "NA"]
        Data.loc[exclude.index,"exclude"] = True 
    
    ## Convert plate columns to phenostring
    Data['pheno_str'] = Data.apply(lambda x:  agar_str(x["agar_AB"],x["agar_A"],x["agar_B"],x["agar_N"]) , axis=1).apply(np.array)    
    Data = Data.drop(columns=['agar_N', 'agar_B', 'agar_A', 'agar_AB'])
    return Data

def transfer_from_to_to_id(row, Data):
    from_id = row["transfer_from_well_id"]
    if type(from_id) == str:
        Data.loc[from_id, "transfer_to_well_id"] = row.name

def r_well_to_id_transfer_infection(Data,strainplate):
    t_end = max(Data.transfer_n)
    data = Data[(Data.transfer_n < t_end)]

    for i, row in data.iterrows():
        t = row["transfer_n"]
        r = row["replicate"]
        
        ## Infection R_well to ID
        ### Only if this well infects someone
        if np.isnan(row.infection_to_well) != True:
            infected_well = Data[(Data["rwell"] == row["infection_to_well"]) & (Data["transfer_n"] == row["transfer_n"]+1) & (Data["replicate"] == row["replicate"]) & (Data["plate"] == row["plate"])]
            Data.loc[i, "infection_to_well_id"] = infected_well.index

        ## Transfer R_well to ID
        ### Turnover
        if type(row.turnover_strain) == str: #WT, A, ... are strings
            pass

        ### No Turnover == Transfer
        else:
            if row["transfer_n"]>0:
                transfer_well = Data[(Data["rwell"] == row["rwell"]) & (Data["transfer_n"] == row["transfer_n"]-1) & (Data["replicate"] == row["replicate"]) & (Data["plate"] == row["plate"])]
                Data.loc[i, "transfer_from_well_id"] = transfer_well.index
        
    return Data

def Infection_Transfer_Turnover_History(Data):
    for t in range(1,max(Data.Time)+1):
        Data_t = Data[Data.Time == t]
        for i, line in Data_t.iterrows():

            ## Check if there was an infection
            if type(line["infected_by_"]) == str: ## Avoid troubles with nans = no infection (nan is a float)
                infecting_id = line["infected_by_"]
                infecting_phenotype = Data.loc[infecting_id, "Phenotype"]
                infecting_contamination_status = Data.loc[infecting_id, "Contaminated"]
                Data.loc[i, "infected_by_"] = infecting_id
                Data.loc[i, "Infecting_Phenotype"] = infecting_phenotype
                Data.loc[i, "Contaminated"] = infecting_contamination_status | Data.loc[i, "Contaminated"]

            ## Check if there was an turnover - Turnovers are already registered. They are only used to reverse identify transfers
            if type(line["Turnover_real"]) == str: ## Avoid troubles with nans = no turnover (nan is a float)
                pass
            ## Otherwise there was an transfer    
            else:
                transfering_id = back_date_id(i) 
                transfering_phenotype = Data.loc[transfering_id, "Phenotype"]
                transfering_contamination_status = Data.loc[transfering_id, "Contaminated"]
                Data.loc[i, "Transfer_Phenotype"] = transfering_phenotype
                Data.loc[i, "Contaminated"] = Data.loc[i, "Contaminated"] | transfering_contamination_status
    return Data


def build_dictionaries():
    ## Build dictionary to assign colors
    hue_infos = {"U":"gray", "S":"green", "A":"blue", "B":"orange", "AB":"red", "A&B":"pink", "Fishy":"black", "Fishy_OD":"Brown", "B_r":"orange", "AB_r":"red", "UI":"gray", "?":"white", "WT":"green", "A_r":"blue"}
    ## Build dictionary to assign strategies
    strategies = {"P1": "No treatment", "P2":"Mono A", "P3":"Mono B","P4":"Combo","P5":"Cycling", "P6":"Mixing"}
    ## Agar Strings to Strains... (AB, A, B, N) -> Binary String -> Decimal -> Strain
    agar_strain_dict = {0:"UI", 5:"A_r", 3:"B_r", 15:"AB_r", 7:"A&B", 1:"WT"}
    ## Well Dictionary
    well = 1
    well_dict = {}
    plate_rows = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]
    for row in plate_rows:
        for col in range(1,25):
            well_dict.update({well : str(row)+str(col)})
            well += 1
    return hue_infos, strategies, plate_rows, agar_strain_dict

def Label_Replicates(Data, plate_rows):
    ## Label Replicates
    for i, line in Data.iterrows():
        row = line["row"]
        row_i = plate_rows.index(row)
        col = line["column"]-1
        replicate = (row_i%2) + col%2*2

        Data.loc[i,"replicate"] = replicate
    return Data

def CreateStrainplates_old_layout_384(T, plate_rows):
    strainplate = pd.DataFrame(columns = ["transfer", "row", "col", "strain_wanted", "strain_real", "contaminated", "replicate"] )
    for t in T:
        for row in plate_rows:
            for col in range(1,25):
                id_ = MakeID(t, "S", row, col)
                if (row == "A") | (row == "B"):
                    if (col == 3) | (col == 4):
                        strainplate.loc[id_, :] = [t, row, col, "A_r", "A_r", False, "?"]
                    elif (col == 5) | (col == 6):
                        strainplate.loc[id_, :] = [t, row, col, "B_r", "B_r", False, "?"]
                    elif (col == 7) | (col == 8):
                        strainplate.loc[id_, :] = [t, row, col, "wt", "wt", False, "?"]
                    elif (col == 9) | (col == 10):
                        strainplate.loc[id_, :] = [t, row, col, "UI", "UI", False, "?"]
                    elif (col == 11) | (col == 12):
                        strainplate.loc[id_, :] = [t, row, col, "bl", "bl", False, "?"]    
                    else:
                        strainplate.loc[id_, :] = [t, row, col, None, None, False, "?"]

                else:
                    strainplate.loc[id_, :] = [t, row, col, None, None, False, "?"]
    return strainplate

def MakeID(t, p, row, col):
    id_ = "t"+str(t)+"_"+str(p)+"_"+str(row)+str(col)
    return id_

def CheckStrainPlates(path, strainplate):
    contaminated_strainplates = pd.read_csv(path+"/pickolo/strainCheckPlates/contaminatedStrainplates.csv")
    for i, entry in contaminated_strainplates.iterrows():
        id_ = MakeID(entry[0], "S", entry[1], entry[2])
        strainplate.loc[id_,"contaminated"] = entry[4]
        strainplate.loc[id_,"strain_real"] = entry[3]
    return strainplate

def get_pathes(experiment):
    # get current directory
    jupyter_path = os.getcwd()

    # prints parent directory
    basepath = os.path.abspath(os.path.join(jupyter_path, os.pardir))

  
    path = basepath + os.sep + "experiments" + os.sep  + experiment 
    return path, os

def get_t_end(instruction_file_names):
    T = []
    for i, row in instruction_file_names.iterrows():
        t = int(row["name"].split("_")[1].split(".")[0])
        T.append(t)
    return max(T)


def AddWellInfos(data, r_well_dict, strategy_dict):
    for i, row in data.iterrows():
        s = i.split("_")
        data.loc[i, "cr_well"] = s[2]
        data.loc[i, "r_well"] = r_well_dict[s[2]]
        data.loc[i, "Strategy"] = strategy_dict[s[1]]
        r = re.sub("\d", "", s[2])
        c = re.sub("\D", "", s[2])
        data.loc[i, "Row"] = r
        data.loc[i, "Column"] = int(c)
        t = re.sub("\D", "", s[0])
        data.loc[i, "Time"] = int(t)
    return data

def Import_Raw_Instructions(path, instruction_file_names):
    
    # read files in the list and append them
    for i, path in instruction_file_names["path"].items(): 
        if i == 0:
            instructions = pd.read_csv(path)
        else:
            instructions = instructions.append(pd.read_csv(path))
    
    # Create the ID column
    t = 't' + instructions["transfer_n"].astype(str)
    p = "_P" + instructions["plate"].astype(str)
    w = "_" + instructions["row"].astype(str) + instructions["col"].astype(str)
    id_ = t + p + w
    
    instructions["ID"] = id_
    instructions = instructions.set_index("ID")

    return instructions


def Label_Replicates(Data, plate_rows):
    ## Label Replicates
    for i, line in Data.iterrows():
        row = line["row"]
        row_i = plate_rows.index(row)
        col = line["col"]-1
        replicate = (row_i%2) + col%2*2

        Data.loc[i,"replicate"] = replicate
    return Data

def include_turnover_from_strainplate(Data, strainplate):

    ## Potential Turnover Strains
    t_strains = ["AB", "A&B", "A_r", "B_r", "WT"]
    Data["Contaminated"] = False

    # Go through each transfer and write the correct wells for each turnover strain to track contamination etc.
    for t in range(0, max(Data.Time)+1):
        for r in [0,1,2,3]:
            Data_t_r = Data[(Data.Time == t) & (Data.Replicate == r)]
            strainplate_t_r = strainplate[(strainplate.Transfer == t) & (strainplate.Replicate == r)]

            # Check where each strain the strainplate is
            for strain in t_strains:
                well =strainplate_t_r[ strainplate_t_r.Strain_wanted == strain]
                indices = Data_t_r[Data_t_r.Turnover == strain].index
                if len(well) != 0: ## If that strain was a turnover strain...
                    Data.loc[indices,"contaminated"] = well["contaminated"].bool()
                    Data.loc[indices,"turnover_well"] = well.index[0]
                    Data.loc[indices,"turnover_real"] = well["strain_real"][0]
    return Data
           

def Include_Instructions_II(Data, instruction_file_names, plate_rows): 
    ## RWell Dictionary
    r_well_dict = {}
    
    for path in instruction_file_names["path"]:
        instructions = pd.read_csv(path)
        for i, line in instructions.iterrows():
            r_well_dict.update({str(line["row"])+str(line["col"]) : line["rwell"]})
        if Data.columns.str.contains("treatment_with").sum(): 
            pass
        else:
            Data.insert(5,"treatment_with", "?")
            Data.insert(8,"infected_by_", "None")
            Data.insert(9,"Infecting_Phenotype" , "None")
            Data.insert(7, "Turnover", "None")

        ## Loop throug instructions for all Plates and figure out who infected whom
        for i, line in  instructions.iterrows():
            p = line["plate"]
            t = line["transfer_n"]
            r = line["row"]
            c = line["col"]
            ## ID that: 
                ## is going to infect another well 
                ## or will be replaced due to turnover 
                ## or is just transfered
            infecting_id = MakeID(t,"P"+str(p),r,c)
                           
            # Document the antibiotic treatment of well
            Data.loc[infecting_id, "treatment_with"] = line["treatment_with"]

            # Document the Turnovers
            if isinstance(line["turnover_strain"], str):
                if line["turnover_strain"] == "wt":
                    Data.loc[infecting_id, "turnover"] = line["turnover_strain"].upper()
                else:
                    Data.loc[infecting_id, "turnover"] = line["turnover_strain"]

            ## Infection
            # Find out which wells belong to infected R-Well
            r_well = line["infection_to_well"]
            Data = Write_Infos_to_Infected_Well(Data, r_well, r_well_dict, infecting_id, plate_rows, t, p)   
    return Data, r_well_dict

def Write_Infos_to_Infected_Well(Data, r_well, r_well_dict, infecting_id, plate_rows, t, p):
    if math.isnan(r_well) is not True and math.isnan(Data.loc[infecting_id, "Replicate"]) == False:
        all_keys = []
        for key, value in r_well_dict.items():
            if(value == r_well):
                 all_keys.append(key)
        ## Find the r_well that belongs to right replicate
        key_dict = {}
        for k in all_keys:
            num = int(k[1:]) 
            r = plate_rows.index(k[0])
            if num%2:
                if r%2:
                    key_dict.update({1:k})
                else:  
                    key_dict.update({0:k})
            else:
                if r%2:
                    key_dict.update({3:k})
                else:
                    key_dict.update({2:k})


        well = key_dict[Data.loc[infecting_id, "Replicate"]]
        ri = re.sub("\d", "", well)
        ci = re.sub("\D", "", well)
        id_infected = MakeID(t,"P"+str(p),ri,ci)

        Data.loc[id_infected, "infected_by_"] = back_date_id(infecting_id) 
        #print("infecting_id", infecting_id, "id_infected", id_infected)
        Data.loc[id_infected, "Infecting_Phenotype"] = Data.loc[infecting_id, "phenotype"]

    return Data

def ImportAgarPlates(file_names, Data, strategies):
    if "agarplates" not in Data.columns:
        Data["agarplates"] = "?"
    for index, file in file_names.iterrows():
        t, p, id_num, ab = DisentangleName(file["name"])
        
        plate_string_dict = {"AB" : 0, "A" : 1, "B" : 2, "N" : 3} 
        
        ## Load single json file: 
        with open(file["path"]) as f:
            agar_plate = json.load(f)
            agar_plate = agar_plate["growth"]
        for well in agar_plate:
            row = well["row"]
            col = well["col"]
            w =  well["well"]
            id_ = MakeID(t, p, row, col)
            Data.loc[id_, "strategy"] = strategies[p]
            ## Manually corrected or auto detected
            
            if len(Data.loc[id_, "agarplates"]) < 4:
                Data.at[id_, "agarplates"] = ["?","?","?","?"]
            plate_string = Data.at[id_, "agarplates"] 
            
            
            if well["manual"] == "NA":
                plate_string[plate_string_dict[ab]] = int(well["auto"])
            else:
                plate_string[plate_string_dict[ab]] = int(well["manual"])
                
            Data.at[id_, "agarplates"]  = plate_string
    return Data

def DisentangleName(name):
    ## 384 Well Plates
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
    
    ## Agar Plates
    else:
        l1 = name.split("_")
        t = re.findall("\d+", l1[0])
        l2 = l1[-1].split(".")
        id_num = l2[0]
        plate_num = l1[1]
        treat = l1[2]
        return t[0], plate_num, id_num, treat
    
    
def back_date_id(id_):
    splitted = id_.split("_")
    t_str = splitted[0]

    t = int(re.sub('\D', '', t_str))

    splitted[0] ="t"+ str(t-1)
    id_old = "_".join(splitted)
    return  id_old

def agar_arr_to_phenotype(agar_arr, agar_strain_dict):
    if len(agar_arr) == 4:
        agar_str = ""
        for a in agar_arr:
            if a:
                agar_str +="1"
            else:
                agar_str += "0"

        agar_num = int(agar_str,2)
        if agar_num in agar_strain_dict.keys():
            strain = agar_strain_dict[agar_num]
        else:
            strain = "Fishy"
    else:
        strain = "?"
    return strain

def EvaluatePhenotype_arr(Data, agar_strain_dict):
    for i in Data.index:
        Data.loc[i, "phenotype"] = agar_arr_to_phenotype(Data.loc[i, "agarplates"], agar_strain_dict)
    return Data

## Count incidence of each Phenotype per Strategy and Replicate
def Summerize_Results_hist(Data):
    summary = Data.groupby(["transfer_n", "replicate", "strategy", "phenotype"]).count()
    summary = summary.reset_index()
    summary = summary.rename(columns={"row": "count", "col" : "fraction"})
    for i, line in summary.iterrows():
        #summary.loc[i, "transfer_n"] = "t" + str(line["transfer_n"]) 
        summary.loc[i, "fraction"] = line["count"]/96*100
    return summary
