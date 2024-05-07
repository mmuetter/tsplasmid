import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.ndimage import label


class KernelPlot:
    def __init__(self, kernel_collector):
        self.parameter_directions = kernel_collector.parameter_directions
        self.kernel_collector = kernel_collector
        self.find_winning_groups()
        self.generate_grid_df()
        self.create_winning_group_dict()
        self.filter_winning_group_dict()
        self.generate_colors()
        self.make_alpha_values()
        self.create_filtered_dataframe()

    def find_winning_groups(self):
        winning_groups = []
        for grid_point in self.kernel_collector.grid_points:
            if not grid_point.empty:
                winning_groups.append(grid_point.winning_group)
        self.winning_groups_all = list(set(winning_groups))

    def generate_colors(self):
        groups = self.winning_group_dict_filtered.keys()
        color_palette = sns.color_palette("Set2", len(groups))
        self.color_dict = dict(zip(groups, color_palette))

    def make_alpha_values(self):
        self.grid_df["alpha"] = round(self.grid_df.num_of_points / 500, 1)
        mask = self.grid_df.alpha > 1
        self.grid_df.loc[mask, "alpha"] = 1
        df = self.grid_df
        x = df['ldaX'].astype(float)
        y = df['ldaY'].astype(float)
        self.x_ticks = np.sort(x.unique())
        self.y_ticks = np.sort(y.unique())
        self.X, self.Y = np.meshgrid(self.x_ticks, self.y_ticks)
        alpha = np.zeros((len(self.y_ticks), len(self.x_ticks)), dtype=float)

        for _, row in df.iterrows():
            ix = np.where(self.x_ticks == row['ldaX'])[0][0]
            iy = np.where(self.y_ticks == row['ldaY'])[0][0]
            alpha[iy, ix] = row['alpha']
        self.alpha_values = alpha

    def generate_grid_df(self):
        self.grid_df = []
        for grid_point in self.kernel_collector.grid_points:
            if not grid_point.empty:
                row = {"ldaX": grid_point.x, "ldaY": grid_point.y,
                       "winning_group": grid_point.winning_group, "num_of_points": grid_point.num_of_points}
                self.grid_df.append(row)
        self.grid_df = pd.DataFrame().from_records(self.grid_df)

    def create_winning_group_dict(self):
        x = self.grid_df['ldaX'].astype(float)
        y = self.grid_df['ldaY'].astype(float)
        self.x_ticks = np.sort(x.unique())
        self.y_ticks = np.sort(y.unique())
        self.X, self.Y = np.meshgrid(self.x_ticks, self.y_ticks)
        self.winning_group_dict = {}
        for group in self.winning_groups_all:
            Z = np.zeros((len(self.y_ticks), len(self.x_ticks)), dtype=bool)
            df = self.grid_df[self.grid_df['winning_group'] == group]
            if len(df) > 3:
                for _, row in df.iterrows():
                    ix = np.where(self.x_ticks == row['ldaX'])[0][0]
                    iy = np.where(self.y_ticks == row['ldaY'])[0][0]
                    Z[iy, ix] = True
                self.winning_group_dict[group] = Z

    def get_sum(self):
        sum_groups = None
        for group in self.winning_group_dict.values():
            if sum_groups is None:
                sum_groups = group.copy()
            else:
                sum_groups = np.logical_or(sum_groups, group)
        self.sum_groups = sum_groups

    def make_island_mask(self):
        self.get_sum()
        labels, _ = label(self.sum_groups)
        big_block_label = np.argmax(np.bincount(labels.flat)[1:]) + 1
        self.island_mask = labels == big_block_label

    def filter_winning_group_dict(self):
        self.make_island_mask()
        self.winning_group_dict_filtered = {}
        for group, array in self.winning_group_dict.items():
            filtered_array = np.logical_and(array, self.island_mask)
            if np.any(filtered_array):
                self.winning_group_dict_filtered[group] = filtered_array

    def contour_plot(self, array, group):
        color = self.color_dict[group]
        array = self.alpha_values * array
        colors_with_alpha = np.zeros((*self.alpha_values.shape, 4))
        colors_with_alpha[:, :, :3] = color
        colors_with_alpha[:, :, 3] = array
        plt.imshow(colors_with_alpha, origin='lower', extent=[
                   self.x_ticks.min(), self.x_ticks.max(), self.y_ticks.min(), self.y_ticks.max()])
        plt.xlabel('ldaX')
        plt.ylabel('ldaY')

    def plot_groups(self, linewidth=0.3):
        self.linewidth = linewidth
        _, self.ax = plt.subplots(figsize=(12, 8))
        for group, array in self.winning_group_dict_filtered.items():
            self.contour_plot(array, group)
        self.legend_handles = [plt.Line2D([], [], marker='o', color='w', markerfacecolor=self.color_dict[label],
                                          markersize=10, label=label) for label in self.winning_group_dict_filtered.keys()]
        self.ax.legend(handles=self.legend_handles, loc='upper center',
                       bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)
        self.ax.set_xlim(self.x_ticks.min(), self.x_ticks.max())
        self.ax.set_ylim(self.y_ticks.min(), self.y_ticks.max())
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(True)
        self.ax.spines['left'].set_visible(True)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')
        self.ax.set_xlabel('ldaX')
        self.ax.set_ylabel('ldaY')
        self.ax.set_title('Kernel Plot')

    def plot_original_axes(self, arrow_length=0.5, color='black', head=0.02):
        for param in self.parameter_directions:
            direction = self.parameter_directions[param].values
            plt.arrow(0, 0, direction[0] * arrow_length, direction[1] * arrow_length,
                      head_width=head, head_length=head, fc=color, ec=color)
            plt.text(direction[0] * arrow_length * 1.2+0.05, direction[1] * arrow_length*1.2,
                     param, color=color, fontsize=15)

    def create_filtered_dataframe(self):
        data = []
        for group, array in self.winning_group_dict_filtered.items():
            indices = np.argwhere(array)
            for index in indices:
                x = self.x_ticks[index[1]]
                y = self.y_ticks[index[0]]
                data.append({'ldaX': x, 'ldaY': y, 'winning_strategy': group})
        self.df = pd.DataFrame(data)
