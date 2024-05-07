import os
import pandas as pd
import shutil
import json
import pickle
import toml
import numpy as np
import matplotlib.colors as mcolors
import colorsys
import seaborn as sns


def random_array(rate, distribution, length):
    array = [None] * length
    for i in range(length):
        if np.random.rand() < rate:
            array[i] = np.random.choice(
                list(distribution.keys()), p=list(distribution.values()))
    return array


def assign_rwell(df, rwell_array, column):
    for i in range(96):
        rwell_idx = df[df["rwell"] == i+1].index
        rwell_value = rwell_array[i]
        df.loc[rwell_idx, column] = rwell_value
    return df


def dict_to_json(path, dictionary):
    with open(path, 'w') as fp:
        json.dump(dictionary, fp, sort_keys=True, indent=4)


def get_base_pathes():
    jupyter_path = os.getcwd()
    basepath = os.path.abspath(os.path.join(jupyter_path, os.pardir))
    return (jupyter_path, basepath)


def get_pathes(experiment):
    jupyter_path = os.getcwd()
    basepath = os.path.abspath(os.path.join(jupyter_path, os.pardir))
    experiment_path = os.path.join(basepath, "experiments", experiment)
    analysis_path = os.path.join(
        basepath, "experiments", experiment, "analysis")
    return jupyter_path, basepath, experiment_path, analysis_path, os


def get_pathes_2(experiment):
    jupyter_path = os.getcwd()
    basepath = os.path.abspath(os.path.join(jupyter_path, os.pardir))
    experiment_path = os.path.join(basepath, "experiments", experiment)
    analysis_path = os.path.join(
        basepath, "experiments", experiment, "analysis")
    obj_path = os.path.join(
        basepath, "experiments", experiment, "analysis", "obj")
    general_obj = os.path.join(basepath, "summary", "obj")
    summary = os.path.join(basepath, "summary")
    tables = os.path.join(basepath, "tables")
    simulations = {"simulations": os.path.join(obj_path, "simulations")}
    return {"jupyter": jupyter_path, "base": basepath, "exp": experiment_path, "analysis": analysis_path, "obj": obj_path, "general_obj": general_obj, "summary": summary, "tables": tables, **simulations}


#  Imports all ... files of a folder and Exports them as a list
def FileNames(ending, path, fold):
    fold_path = path+"/"+fold
    tmp = os.listdir(fold_path)
    files = pd.DataFrame(columns=["name", "path"])
    for file in tmp:
        if ending in file:
            files.loc[len(files), "name"] = file
            files.loc[len(files)-1, "path"] = fold_path+"/"+file
    # Exclude hidden lock files
    files = files[~files.name.str.contains("lock")]
    return files


def get_t_end(instruction_file_names):
    T = []
    for i, row in instruction_file_names.iterrows():
        t = int(row["name"].split("_")[1].split(".")[0])
        T.append(t)
    return max(T)


def multiply_setup_plates(path, os):
    agar_file_names = FileNames(".json", path, "pickolo"+os.sep+"json")
    setup_plate_names = []
    for name in agar_file_names["name"]:
        if name[0:2] == "t0":
            setup_plate_names.append(name)
    setup_plate_names

    ending = ".json"
    folder = "pickolo"+os.sep+"json"
    json_path = path+os.sep + folder

    for name in setup_plate_names:
        if "setup" not in name:
            split_name = name.split("_")
            for i in range(2, 7):
                p = "P" + str(i)
                new_split = split_name
                new_split[1] = p
                new_split[3] = "setup.json"
                new_name = "_".join(new_split)
                src = json_path + os.sep + name
                target = json_path + os.sep + new_name
                shutil.copyfile(src, target)


def obj_to_pickle(name, data):
    with open(name, 'wb') as f:
        pickle.dump(data, f)


def load_pickle_obj(name):
    with open(name, 'rb') as f:
        data = pickle.load(f)
    return data


def save_for_latex(df, math_cols, path, round=2):
    df = df.reset_index()
    columns = set(df.columns)
    c = columns.difference(set(math_cols))
    for col in math_cols:
        df[col] = df[col].apply(
            lambda x: adapt_entry(x, make_math=True, r=round))
    for col in c:
        df[col] = df[col].apply(lambda x: adapt_entry(x, r=round))
    df.columns = [w.replace('_', '') for w in df.columns]
    df.columns = [w.replace('&', 'u') for w in df.columns]
    df.to_csv(path, sep=";")
    print(df)


def adapt_entry(x, r=2, make_math=False):
    if type(x) == float:
        x = round(x, r)
    elif make_math:
        x = "$"+x+"$"
    elif x == "A&B":
        x = x.replace("A&B", "(A\\_r\\&B\\_r)")
    elif type(x) == str:
        x = x.replace("_", "\\_")
    elif (type(x) == int) | (type(x) == bool):
        pass
    else:
        print("wtf?", x)
    return x


def load_json(path):
    f = open(path)
    file = json.load(f)
    f.close()
    return file


def dump_toml(name, data):
    f = open(name, 'w')
    toml.dump(data, f)
    f.close()


def load_toml(name):
    return toml.load(name)


def build_dictionaries_2():
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

    rename = {
        "wt": "S",
        "WT": "S",
        "UI": "U",
        "Fishy": "Other"
    }
    return {"strategies": strategies, "rows": plate_rows, "agar_strains": agar_strain_dict, "rename": rename}


def build_dictionaries_3():
    # Build dictionary to assign colors
    hue_infos = {"U": "gray", "S": "green", "A": "blue", "B": "orange", "AB": "red", "A&B": "pink", "Fishy": "black",
                 "Fishy_OD": "Brown", "B_r": "orange", "AB_r": "red", "UI": "gray", "?": "white", "WT": "green", "A_r": "blue"}
    #  Build dictionary to assign strategies
    strategies = {"P1": "No treatment", "P2": "Mono A",
                  "P3": "Mono B", "P4": "Combination", "P5": "Cycling", "P6": "Mixing"}
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

    rename = {
        "wt": "S",
        "WT": "S",
        "UI": "U",
        "Fishy": "Other"
    }

    axs_dict = {"No treatment": (0, 0), "Mono A": (0, 1), "Mono B": (
        0, 2), "Mixing": (1, 0), "Cycling": (1, 1), "Combination": (1, 2)}
    return {"plates": strategies, "rows": plate_rows, "agar_strains": agar_strain_dict, "rename": rename, "axs_dict": axs_dict}


def load_experiment(exp, src="01b_Data"):
    pathes = get_pathes_2(exp+"_output")
    dictionaries = build_dictionaries_2()
    Data = pd.read_pickle(pathes["obj"]+os.sep + src + ".pkl")
    Data["exp"] = exp
    t_end = max(Data["transfer_n"])+1
    with open(os.path.join(pathes["exp"], "exp_pars.json")) as json_data:
        exp_pars = json.load(json_data)
    strains = ["U", "S", "A_r", "B_r", "A&B", "AB_r"]
    return strains, exp_pars, dictionaries, pathes, Data, t_end


def rename_strings(data, rename_dict):
    def replace_string(x):
        if isinstance(x, str):
            if x in rename_dict:
                return rename_dict[x]
        elif isinstance(x, (list, tuple, np.ndarray)):
            return [replace_string(item) for item in x]
        return x

    return data.applymap(replace_string)


def removekey(d, key):
    r = dict(d)
    del r[key]
    return r


def sort_dictionary(dictionary, order):
    c = {}
    for _, o in enumerate(order):
        if o in dictionary.keys():
            c.update({o: dictionary[o]/94})
        else:
            c.update({o: 0})
    return c


def df_to_dict(df, i_col=0):
    return dict(zip(df.index, df.iloc[:, i_col].values))


def protocol_anova(phenos, anovas):
    anova_results = []
    for pheno, anova in zip(phenos, anovas):
        label = pheno.split("_")
        label = " ".join(label)
        print(label)
        anova_results.append(
            {"group": label, "statistic": anova.statistic, "Pvalue": anova.pvalue, "significant": anova.pvalue < 0.05})
    anova_results = pd.DataFrame().from_records(anova_results)
    return anova_results


def save_tukey_results(tuk, exp, label, path):
    summary = tuk.summary()
    name = "tukey_"+exp+"_"+label+".csv"
    file_path = os.path.join(path, name)
    data = summary.data
    headers = data[0]
    df = pd.DataFrame(data[1:], columns=headers)
    df.columns = ["groupI", "groupII", "meandiff",
                  "pAdj", "lower", "upper", "reject"]
    df["pValue"] = tuk.pvalues
    print(file_path)
    save_for_latex(df, [], file_path, round=4)


def adjust_color_for_contrast(color, saturation_factor=1.5, contrast = 1, background_color='white'):
    # Convert color to HSV (Hue, Saturation, Value)
    hsv = colorsys.rgb_to_hsv(*mcolors.to_rgb(color))
    
    # Increase saturation
    new_saturation = min(hsv[1] * saturation_factor, 1)
    
    # Adjust value based on background color
    if background_color == 'white':
        
        new_value = min(hsv[2] /contrast, 1)  # Darken for white background
    else:
        new_value = min(hsv[2] * contrast, 1)  # Brighten for black background

    # Convert back to RGB
    new_color = colorsys.hsv_to_rgb(hsv[0], new_saturation, new_value)
    return mcolors.to_hex(new_color)

   