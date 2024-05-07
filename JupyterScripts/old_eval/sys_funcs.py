import numpy as np
import pandas as pd
import pickle

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as file:
        pickle.dump(obj, file, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as file:
        return pickle.load(file)