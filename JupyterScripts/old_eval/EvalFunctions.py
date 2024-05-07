import os
import pandas as pd
import re
import json
import xml.etree.ElementTree as ET
import re
import numpy as np
import math


### Imports all ... files of a folder and Exports them as a list
def FileNames(ending, path, fold):
    fold_path = path+"/"+fold
    tmp = os.listdir(fold_path)
    files = pd.DataFrame(columns = ["name", "path"])
    for file in tmp:
        if ending in file:
            files.loc[len(files),"name"] = file
            files.loc[len(files)-1,"path"] = fold_path+"/"+file  
    #Exclude hidden lock files        
    files = files[~files.name.str.contains("lock")]
    return files

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

def MakeID(t, p, row, col):
    id_ = "t"+str(t)+"_"+str(p)+"_"+str(row)+str(col)
    return id_

def ImportAgarPlates(file_names, Data, strategies):
    for index, file in file_names.iterrows():
        t, p, id_num, ab = DisentangleName(file["name"])

        ## Load single json file: 
        with open(file["path"]) as f:
            agar_plate = json.load(f)
            agar_plate = agar_plate["growth"]
        for well in agar_plate:
            row = well["row"]
            col = well["col"]
            w =  well["well"]
            id_ = MakeID(t, p, row, col)
            Data.loc[id_, "Well"] = w 
            Data.loc[id_, "Row"] = well["row"]
            Data.loc[id_, "Time"] = t
            Data.loc[id_, "Column"] = col
            Data.loc[id_, "Strategy"] = strategies[p]
            if well["manual"] == "NA":
                Data.loc[id_, "Agar_"+ab] = well["auto"]
            else:
                Data.loc[id_, "Agar_"+ab] = well["manual"]
    Data["Time"] = pd.to_numeric(Data["Time"])
    return Data

def EvaluatePhenotype(Data):
    for i, data in Data.iterrows():
        if data["Agar_N"] == 0 and data["Agar_A"] + data["Agar_B"] + data["Agar_AB"] > 0:
            Data.loc[i,"Phenotype"] = "Fishy"
        elif data["Agar_AB"] == 1 and data["Agar_N"] + data["Agar_A"]+ data["Agar_B"] != 3:
            Data.loc[i,"Phenotype"] = "Fishy"
        elif data["Agar_AB"] == 1 :
            Data.loc[i,"Phenotype"] = "AB"
        elif data["Agar_A"] == 1 and data["Agar_B"] == 1:
            Data.loc[i,"Phenotype"] = "A&B"
        elif data["Agar_A"] == 1 :
            Data.loc[i,"Phenotype"] = "A"
        elif data["Agar_B"] == 1 :
            Data.loc[i,"Phenotype"] = "B"
        elif data["Agar_N"] == 1 :
            Data.loc[i,"Phenotype"] = "S"
        elif data["Agar_N"] + data["Agar_A"] + data["Agar_B"] +data["Agar_AB"] == 0 :
            Data.loc[i,"Phenotype"] = "U"
        else:
            Data.loc[i,"Phenotype"] = "Fishy"
    return Data

def LoadODs(od_file_names, Data):
    for name, path in zip(od_file_names["name"], od_file_names["path"]):

        # Scan XML for ODs
        tree = ET.parse(path)
        root = tree.getroot()
        it = root.iter('Well')
        t, barcode, p, m = DisentangleName(name)
        
        ## First Plate (t0 B) are not read. Correct Data using the min OD Value (at this time nothing is grown either way)
        ind = Data.index.str.contains("t0")
        min_OD = Data["OD_A"].min()
        Data.loc[ind, "OD_B"] = min_OD
        
        for i in it:
            position = i.attrib["Pos"]
            temp = re.compile("([a-zA-Z]+)([0-9]+)")
            res = temp.match(position).groups()
            row = res[0]
            col = res[1]
            for od in i:
                id_ = MakeID(t, p, row, col)
                Data.loc[id_, "Time"] = t


                ## A... After, B ... Before : Same Well, Same Plate: t0 (B) -> t1 (A): correct to t0
                if m == "AI" :
                    Data.loc[id_, "OD_A"] = float(od.text)
                elif m == "BI":                                        
                    Data.loc[id_, "OD_B"] = float(od.text)
                else:
                    print("Something is wrong...")
            Data.loc[id_, "Row"] = row
            Data.loc[id_, "Col"] = col

        
    return Data

def CalculateGrowthAndCleanDataframe(Data):
    val = np.nan
    for i, row in Data.iterrows(): 
        # Remove incomplete Rows (e.g. No second OD measure or no Agar Plate Measures)
        if row.isnull().sum() > 1:
            Data = Data.drop(i)

        # If row complete: Calculate if growth in well happend
        else:
            rel = (row["OD_A"]-row["OD_B"])/row["OD_B"]
            if rel > 2:
                Data.loc[i, "Growth"] =  True
            else:
                Data.loc[i, "Growth"] = False
    return Data

def Label_Replicates(Data, plate_rows):
    ## Label Replicates
    for i, line in Data.iterrows():
        row = line["Row"]
        row_i = plate_rows.index(row)
        col = line["Column"]-1
        replicate = (row_i%2) + col%2*2

        Data.loc[i,"Replicate"] = replicate
    return Data

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
        Data.loc[id_infected, "Infecting_Phenotype"] = Data.loc[infecting_id, "Phenotype"]

    return Data

def Include_Instructions(Data, instruction_file_names, plate_rows):
    Data["Data_complete"] = False
   
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
            Data.insert(10,"infected_by_", "None")
            Data.insert(11,"Infecting_Phenotype" , "None")
            Data.insert(12, "Turnover", "None")

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
                Data.loc[infecting_id, "Turnover"] = line["turnover_strain"]

            ## Infection
            # Find out which wells belong to infected R-Well
            r_well = line["infection_to_well"]
            Data = Write_Infos_to_Infected_Well(Data, r_well, r_well_dict, infecting_id, plate_rows, t, p)   
    Data["Data_complete"] = Data["Agar_A"].isnull()
    return Data, r_well_dict

## Count incidence of each Phenotype per Strategy and Replicate
def Summerize_Results_hist(Data):
    summary = Data.groupby(["Time", "Replicate", "Strategy", "Phenotype"]).count()
    summary = summary.reset_index()
    summary = summary.rename(columns={"Row": "Count", "Column" : "Fraction"})
    for i, line in summary.iterrows():
        summary.loc[i, "Time"] = "t" + str(line["Time"]) 
        summary.loc[i, "Fraction"] = line["Count"]/96*100
    return summary


def back_date_id(id_):
    splitted = id_.split("_")
    t_str = splitted[0]

    t = int(re.sub('\D', '', t_str))

    splitted[0] ="t"+ str(t-1)
    id_old = "_".join(splitted)
    return  id_old
