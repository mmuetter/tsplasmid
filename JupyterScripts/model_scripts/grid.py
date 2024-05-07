import pandas as pd
import numpy as np
from scipy.spatial import distance
import pandas as pd
import os
import pickle


class GridPoint:
    def __init__(self, x, y, d50, data, strategy_columns):
        self.x = x
        self.y = y
        self.d50 = d50
        self.data = data.copy()
        self.strategy_columns = strategy_columns[:5]
        self.get_distance()
        self.keep_closest_points()
        self.get_weights()
        self.make_weighted_averages()

    def get_distance(self):
        lda_axis_1 = self.data["LDA_axis_1"]
        lda_axis_2 = self.data["LDA_axis_2"]
        distances = distance.cdist(
            [[self.x, self.y]], np.column_stack((lda_axis_1, lda_axis_2)))
        self.data["distance"] = distances[0, 0:]

    def keep_closest_points(self):
        threshold = self.data['distance'].quantile(0.01)
        mask = self.data['distance'] <= threshold
        self.data = self.data[mask]

    def get_weights(self):
        self.data["weights"] = self.data.apply(
            lambda x: self.weight(x.distance), axis=1)
        self.sum_of_weights = self.data.weights.sum()

    def weight(self, d):
        a = -np.log(0.5)/(self.d50**2)
        weight = np.exp(-a*(d**2))
        return weight

    def make_weighted_averages(self):
        self.weighted_averages = {}
        best_f = 0
        for col in self.strategy_columns:
            weighted_av = (
                self.data[col] * self.data.weights).sum()/self.sum_of_weights
            self.weighted_averages.update(
                {col: weighted_av})
            if best_f < weighted_av:
                best_f = weighted_av
                self.winning_strategy = col


class MakeGrid:
    def __init__(self, sensitivity, path, d50=0.025):
        self.d50 = d50
        self.strategy_columns = sensitivity.strategies
        self.parameter_directions = sensitivity.parameter_directions
        self.path = path
        self.file_path = os.path.join(
            path, "grid" + sensitivity.name_add + ".pkl")

    def create_grid_points(self, sensitivity, num_ticks):
        self.num_ticks = num_ticks
        lda_axis_1 = sensitivity.evaluations_df["LDA_axis_1"]
        lda_axis_2 = sensitivity.evaluations_df["LDA_axis_2"]
        min_x, max_x = min(lda_axis_1), max(lda_axis_1)
        min_y, max_y = min(lda_axis_2), max(lda_axis_2)
        self.x_vals = np.linspace(min_x, max_x, self.num_ticks)
        self.y_vals = np.linspace(min_y, max_y, self.num_ticks)
        self.grid_points = []
        for x in self.x_vals:
            print(x)
            for y in self.y_vals:
                grid_point = GridPoint(
                    x, y, self.d50, sensitivity.evaluations_df, self.strategy_columns)
                self.grid_points.append(grid_point)

    def summarize(self):
        df = []
        for grid_point in self.grid_points:
            row = {"sum_of_weights": grid_point.sum_of_weights,
                   "winning_strategy": grid_point.winning_strategy,
                   "ldaX": grid_point.x,
                   "ldaY": grid_point.y}
            df.append(row)
        self.summary = pd.DataFrame().from_records(df)

    def make_alpha(self,  upper=1, lower=0.1):
        df = self.summary.copy()
        df["alpha"] = None
        df.alpha = round((df.sum_of_weights-lower)/(upper-lower), 2)
        df.loc[df.alpha < 0, "alpha"] = 0
        df.loc[df.alpha > 1, "alpha"] = 1
        self.summary["alpha"] = df.alpha

    def save(self):
        with open(self.file_path, 'wb') as file:
            pickle.dump(self, file)

    def load(self):
        with open(self.file_path, 'rb') as file:
            obj = pd.read_pickle(file)
        self.__dict__.update(obj.__dict__)
