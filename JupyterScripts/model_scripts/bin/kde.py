import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D
import seaborn as sns


class kdeObj:
    def __init__(self, sensitivity_analysis, num_points=1000):
        self.sensitivity_analysis = sensitivity_analysis
        self.lda_df = sensitivity_analysis.lda_df
        self.parameters = sensitivity_analysis.parameter_names
        self.parameter_directions = sensitivity_analysis.parameter_directions
        self.groups = sensitivity_analysis.groups
        self.num_points = num_points
        self.get_grid()
        self.get_kde_data()
        colors = sns.color_palette('tab10', len(self.groups) + 1)
        self.strategy_colors = dict(zip(self.groups, colors))

    def get_grid(self):
        x_min, x_max = self.lda_df['LDA Axis 1'].min(
        ), self.lda_df['LDA Axis 1'].max()
        y_min, y_max = self.lda_df['LDA Axis 2'].min(
        ), self.lda_df['LDA Axis 2'].max()
        x_grid, y_grid = np.mgrid[x_min:x_max:self.num_points *
                                  1j, y_min:y_max:self.num_points*1j]
        self.grid_points = np.vstack([x_grid.ravel(), y_grid.ravel()])

    def get_kde(self, strategies, bw_method=.3):
        df = self.lda_df
        df = df[df.winning_group.isin(strategies)]
        if len(df) > 2:
            x, y = df['LDA Axis 1'], df['LDA Axis 2']
            kde = gaussian_kde([x, y], bw_method=bw_method)
            return kde.evaluate(self.grid_points) * len(df)
        else:
            return []

    def get_frequencies(self, density_thresh=100):
        total = 0 * self.kde_data["Combo"]
        self.frequencies = {}
        for strategy, kde in self.kde_data.items():
            total += kde
        self.point_density = total
        self.density_mask = total > density_thresh
        for strategy, kde in self.kde_data.items():
            frequencies = np.where(total >= density_thresh, kde / total, 0)
            self.frequencies.update({strategy: frequencies})

    def get_kde_data(self):
        self.kde_data = {}
        for group in self.groups:
            kde = self.get_kde([group])
            if len(kde) > 0:
                self.kde_data.update({group: kde})

    def plot_kde(self, strategy, ax=plt):
        self.ax = ax
        alpha = self.kde_data[strategy]/len(self.kde_data[strategy])
        self.ax = plt.scatter(
            self.grid_points[0], self.grid_points[1], alpha=alpha, color=self.strategy_colors[strategy])

    def plot_frequency(self, strategy, thresh_low, thresh_high):
        frequencies = self.frequencies[strategy].copy()
        color = self.strategy_colors[strategy]

        fill_contours = None

        # Plot filled contour if frequencies exceed thresh_high
        if np.any(frequencies > thresh_high):
            fill_contours = plt.contourf(self.grid_points[0].reshape((self.num_points, self.num_points)),
                                         self.grid_points[1].reshape(
                                             (self.num_points, self.num_points)),
                                         frequencies.reshape(
                                             (self.num_points, self.num_points)),
                                         levels=[thresh_high, 1.0],
                                         colors=[(*color[:3], 0.5)])

        # Plot contour if frequencies exceed thresh_low
        if np.any(frequencies > thresh_low):
            contour = plt.contour(self.grid_points[0].reshape((self.num_points, self.num_points)),
                                  self.grid_points[1].reshape(
                                      (self.num_points, self.num_points)),
                                  frequencies.reshape(
                                      (self.num_points, self.num_points)),
                                  levels=[thresh_low],
                                  colors=[(*color[:3], 1.0)])
            if fill_contours is None:
                fill_contours = contour

        return fill_contours

    def plot_frequencies(self, thresh_low=0.25, thresh_high=0.25):
        self.fig, self.ax = fig, ax = plt.subplots(figsize=(12, 12))
        legend_elements = []
        for group in self.frequencies.keys():
            self.plot_frequency(group, thresh_low, thresh_high)
            legend_elements.append(
                Line2D([0], [0], color=self.strategy_colors[group], lw=2, label=group))
        ax.legend(handles=legend_elements, title='Best Strategy')
        return fig, ax

    def plot_original_axes(self, arrow_length=2.5, color='black'):
        for param in self.parameters:
            direction = self.parameter_directions[param].values
            self.ax.arrow(0, 0, direction[0] * arrow_length, direction[1] * arrow_length,
                          head_width=0.1, head_length=0.1, fc=color, ec=color)
            self.ax.text(direction[0] * arrow_length * 1.3, direction[1]
                         * arrow_length * 1.3, param, color=color, fontsize=15)

    def transform_point(self, point_params):
        lda_directions = self.parameter_directions
        transformed_point = np.zeros(2)

        for param, value in point_params.items():
            direction = lda_directions[param].values
            transformed_point += value * direction
        return transformed_point

    def plot_point(self, point_params, marker='x', color='red', label='Transformed Point'):
        point = self.transform_point(point_params)
        plt.scatter(point[0], point[1], color=color,
                    marker=marker, s=100, label=label)
        plt.annotate(label, (point[0], point[1]), xytext=(
            10, 10), textcoords='offset points', fontsize=12)
