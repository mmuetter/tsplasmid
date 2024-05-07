import pandas as pd
import numpy as np
import os
import re
from datetime import datetime


#  Excluded Plates (Dropped etc)
# Put exclude_plates.csv in the Notes folder with columns:[plate, transfer] to define what plates should be excluded
def exclude_unusable_plates(Data, path):
    exclude = pd.read_csv(path+os.sep+"Notes"+os.sep+"exclude_plates.csv")
    P_ex = list(exclude["P_ex"])
    t_ex = list(exclude["t_ex"])
    Data["exclude"] = False
    Data["contaminated"] = False

    for p, t in zip(P_ex, t_ex):
        print(p, t)
        Data.loc[(Data.plate == p) & (Data.transfer_n == t), "exclude"] = True
        Data.loc[(Data.plate == p) & (Data.transfer_n == t),
                 "comment"] = "See notes"
    return Data


def added_wells(x):
    w = []
    if type(x.turnover_strain) == str:
        w.append(x.turnover_id)
    if bool(x.received_transfer_strain):
        w.append(x.transfer_from_well_id)
    if bool(x.infected_by_wells):
        w.extend(x.infected_by_wells)
    return w


def added_strains(x):
    w = []
    if type(x.turnover_strain) == str:
        w.append(x.turnover_strain)
    if bool(x.received_transfer_strain):
        w.append(x.received_transfer_strain)
    if bool(x.infected_by_strains):
        w.extend(x.infected_by_strains)
    return w


def add_time_stamps(Data, path):
    t_path = os.path.join(path, "times.csv")
    time_table = pd.read_csv(t_path)
    barcodes = time_table[time_table["barcode"].str.isnumeric()
                          ]["barcode"].unique()
    events = ['turnover_start', 'incubation_start',
              'replication_start', 'incubation_end']
    relativ_times = ['t_turnover_start', 't_incubation_start',
                     't_replication_start', 't_incubation_end']

    #event_times = pd.DataFrame(columns=["barcode", "transfer", "plate"]+events)
    event_times = []
    for bc in barcodes:
        bc_plate = time_table[time_table.barcode == bc]
        row = {"plate": bc_plate["plate"].unique(
        )[0], "barcode": bc, "transfer": min(bc_plate.transfer)}
        for event in events:
            time = bc_plate.loc[bc_plate.event == event, "time"].to_numpy()
            l, = time.shape
            if l > 0:
                time_array = np.array(
                    re.split('-|:|T|Z', time[0]))[0:-1].astype(int)
                d_time = datetime(year=time_array[0], month=time_array[1], day=time_array[2],
                                  hour=time_array[3], minute=time_array[4], second=time_array[5])
                row.update({event: d_time})
            else:
                row.update({event: None})

        event_times.append(row)
    event_times = pd.DataFrame().from_records(event_times)

    t0 = min(list(event_times.turnover_start.dropna(
    ))+list(event_times.incubation_start.dropna())+list(event_times.replication_start.dropna()))
    for rel_t, event in zip(relativ_times, events):
        tmp = event_times[event].dropna()
        event_times.loc[tmp.index, rel_t] = event_times.loc[tmp.index, :].apply(
            lambda x: (x[event] - t0).days*24 + (x[event] - t0).seconds/3600, axis=1)
    event_times.head()
    time_cols = events+relativ_times+["barcode"]

    for i, row in event_times.iterrows():
        Slice = Data[(Data.plate == row["plate"]) & (
            Data.transfer_n == row["transfer"])]
        for col in time_cols:
            Data.loc[Slice.index, col] = row[col]
    return Data
