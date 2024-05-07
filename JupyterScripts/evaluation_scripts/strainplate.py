import pandas as pd


def CreateStrainplate_layout_384(T, plate_rows):
    strainplate = pd.DataFrame(columns=[
                               "transfer", "row", "col", "strain_wanted", "strain_real", "contaminated", "replicate"])
    for t in T:
        for row in plate_rows:
            for col in range(1, 25):
                id_ = "t"+str(t)+"_S_"+row+str(col)
                if (row == "G") | (row == "H"):
                    if (col == 9) | (col == 10):
                        strainplate.loc[id_, :] = [
                            t, row, col, "A_r", "A_r", False, "?"]
                    elif (col == 14) | (col == 13):
                        strainplate.loc[id_, :] = [
                            t, row, col, "B_r", "B_r", False, "?"]
                    elif (col == 5) | (col == 6):
                        strainplate.loc[id_, :] = [
                            t, row, col, "wt", "wt", False, "?"]
                    elif (col == 17) | (col == 18):
                        strainplate.loc[id_, :] = [
                            t, row, col, "AB_r", "AB_r", False, "?"]
                    elif (col == 21) | (col == 22):
                        strainplate.loc[id_, :] = [
                            t, row, col, "UI", "UI", False, "?"]
                    elif (col == 23) | (col == 24):
                        strainplate.loc[id_, :] = [
                            t, row, col, "bl", "bl", False, "?"]
                    else:
                        strainplate.loc[id_, :] = [
                            t, row, col, None, None, False, "?"]

                else:
                    strainplate.loc[id_, :] = [
                        t, row, col, None, None, False, "?"]
    return strainplate
