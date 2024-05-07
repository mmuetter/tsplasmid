import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, f_oneway
import seaborn as sns
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def norm_second_best(row, strategies):
    for strategy in strategies:
        row[strategy] = (row[strategy] - row['second_best_value']
                         ) / row['second_best_value']
    return row


class PrepareSeperation:
    def __init__(self, df, strategies=["Combo", "Cycling", "Mixing", "Mono A", "Mono B", "No treatment"], parameters=["turnover", "infection", "cA", "cB", "cS", "cU"]):
        self.df = df
        self.strategies = strategies
        self.parameters = parameters
        self.columns = self.strategies + self.parameters
        self.df = self.df[self.columns]
        self.get_relative_quality()
        self.result = self.quality_df[self.quality_df["winner"] > 0.05]

    def get_relative_quality(self, method=norm_second_best):
        self.quality_df = self.df.copy()
        self.quality_df['second_best'] = self.quality_df[self.strategies].apply(
            lambda row: row.nlargest(2).index[-1], axis=1)
        self.quality_df['second_best_value'] = self.quality_df.apply(
            lambda row: row[row['second_best']], axis=1)
        self.quality_df['best_strategy'] = self.quality_df[self.strategies].idxmax(
            axis=1)
        self.quality_df = self.quality_df.apply(
            lambda x: method(x, self.strategies), axis=1)
        self.quality_df['winner'] = self.quality_df[self.strategies].max(
            axis=1)


class Simulation:
    def __init__(self, df, par):
        self.par = par
        self.sim_df = df
        self.make_tEnd()
        self.get_fU()
        self.add_parameter()
        self.results = self.results.set_index("rep")

    def add_parameter(self):
        for rep, parameter_set in enumerate(self.par):
            df = self.results
            mask = df.rep == rep
            df.loc[mask, "turnover"] = parameter_set["turnover"]
            df.loc[mask, "infection"] = parameter_set["infection"]
            df.loc[mask, "cA"] = parameter_set["community"]["A_r"]
            df.loc[mask, "cB"] = parameter_set["community"]["B_r"]
            df.loc[mask, "cS"] = parameter_set["community"]["S"]
            df.loc[mask, "cU"] = parameter_set["community"]["U"]
            self.results = df

    def make_tEnd(self):
        self.t_end = t_end = self.sim_df.transfer_n.max()
        self.sim_df_end = self.sim_df[self.sim_df.transfer_n.isin(
            [t_end, t_end-1, t_end-2, t_end-3])].copy()
        self.sim_df_end_U = self.sim_df_end[self.sim_df_end.x_hat.isin(["U"])]
        columns = ["strategy",  "rep", "f"]
        self.sim_df_end_U = self.sim_df_end_U[columns].groupby(
            columns[:-1]).mean().reset_index()

    def get_fU(self):
        df = self.sim_df_end_U
        pivoted_df = df.pivot(index='rep', columns='strategy', values='f')
        pivoted_df = pivoted_df.reset_index()
        self.results = pivoted_df


class StrategyComparison:
    def __init__(self, simulation: pd.DataFrame, strategy1: str, strategy2: str):
        self.strategy1, self.strategy2 = strategy1, strategy2
        columns = [strategy1, strategy2, "turnover", "infection", "cA", "cB"]
        self.data = simulation[columns].copy()
        self.data["difference"] = simulation[strategy1] - simulation[strategy2]
        self.conduct_anova()
        self.conduct_pearson()

    def conduct_anova(self, alpha=0.05):
        df = self.data
        grouped_data = [df[df['infection'] == rate]['difference']
                        for rate in df['infection'].unique()]
        f_statistic, p_value = f_oneway(*grouped_data)
        self.anova = {"f_statistic": f_statistic, "p_value": p_value}
        self.anova_significant = p_value < alpha

    def conduct_pearson(self, alpha=0.05):
        corr_coef, p_value = pearsonr(
            self.data['infection'], self.data['difference'])
        self.pearson = {"corr_coef": corr_coef, "p_value": p_value}
        self.pearson_significant = p_value < alpha


def create_heatmap(df, x_col, y_col, color_col):
    x_values = df[x_col]
    y_values = df[y_col]
    color_values = df[color_col]
    unique_x = np.unique(x_values)
    unique_y = np.unique(y_values)
    heatmap = np.zeros((len(unique_y), len(unique_x)))

    for r, x_val in enumerate(unique_x):
        for c, y_val in enumerate(unique_y):
            idx = (x_values == x_val) & (y_values == y_val)
            if any(idx):
                heatmap[c, r] = color_values[idx].iloc[0]

    return heatmap, unique_x, unique_y


def plot_heatmap(df, heatmap, unique_x, unique_y, color_col, x_col, y_col):
    _, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    im = ax.imshow(heatmap, cmap='RdYlGn', aspect='auto')

    cbar = ax.figure.colorbar(im, ax=ax, pad=0.01)

    color_col_type = df[color_col].dtype
    if color_col_type == bool:
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['False', 'True'])
        cbar.set_label(color_col)
    elif np.issubdtype(color_col_type, np.number):
        cbar.set_label(color_col)

    ax.set_xticks(np.arange(len(unique_x)))
    ax.set_yticks(np.arange(len(unique_y)))
    ax.set_xticklabels(unique_x)
    ax.set_yticklabels(unique_y)

    ax.set_xticks(np.arange(len(unique_x)) + 0.5, minor=True)
    ax.set_yticks(np.arange(len(unique_y)) + 0.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"Heatmap of {x_col} and {y_col} with {color_col}")

    plt.show()


def plot_heatmap_from_dataframe(df, x_col, y_col, color_col):
    heatmap, unique_x, unique_y = create_heatmap(df, x_col, y_col, color_col)
    plot_heatmap(df, heatmap, unique_x, unique_y, color_col, x_col, y_col)


def make_comparison(pars, par_space, simulations, strategy1, strategy2):
    results = []
    comparisons = []
    heatmap_df = []
    sims = []

    for par, pair, sim_df in zip(pars, par_space, simulations):
        sim = Simulation(sim_df, par)
        results.append(sim.results)
        comparison = StrategyComparison(sim.results, strategy1, strategy2)
        comparisons.append(comparison)
        pair.update({"significance anova": comparison.anova_significant,
                     "significance pearson": comparison.pearson_significant,
                     "pearson pvalue": comparison.pearson["p_value"],
                     "anova pvalue": comparison.anova["p_value"],
                     "pearson correlation coef": comparison.pearson["corr_coef"],
                     "comparison": comparison})
        heatmap_df.append(pair)
        sims.append(sim)
    heatmap_df = pd.DataFrame().from_records(heatmap_df)

    return heatmap_df, comparisons, sims


class LDA:
    def __init__(self, prepare_obj):
        self.quality_df = prepare_obj.quality_df
        self.strategies = prepare_obj.strategies
        colors = sns.color_palette('pastel', len(self.strategies))
        self.strategy_colors = dict(zip(self.strategies, colors))

        self.parameters = prepare_obj.parameters

        self.compute_lda()

    def compute_lda(self):
        # Fit and transform the original data using LDA
        self.lda = lda = LinearDiscriminantAnalysis(n_components=2)
        lda_data = lda.fit_transform(
            self.quality_df[self.parameters], self.quality_df['best_strategy'])

        # Compute correlations between LDA axes and original parameter columns
        correlations = np.corrcoef(
            lda_data.T, self.quality_df[self.parameters].T)

        self.lda_df = pd.DataFrame(
            lda_data, columns=['LDA Axis 1', 'LDA Axis 2'])
        self.lda_df['best_strategy'] = self.quality_df['best_strategy']
        self.parameter_directions = pd.DataFrame(
            correlations[:2, 2:], columns=self.parameters)

    def plot_lda(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        for strategy in self.strategies:
            strategy_data = self.lda_df[self.lda_df['best_strategy'] == strategy]
            ax.scatter(strategy_data['LDA Axis 1'], strategy_data['LDA Axis 2'],
                       color=self.strategy_colors[strategy], label=strategy)
        ax.set_xlabel('LDA Axis 1')
        ax.set_ylabel('LDA Axis 2')
        ax.set_title('Preferable Strategy Map')
        ax.legend(title='Best Strategy')
        self.ax = ax
        self.fig = fig
        self.plot_original_axes()

    def plot_original_axes(self, arrow_length=2.5):
        for param in self.parameters:
            color = 'black'  # Set the color to black
            direction = self.parameter_directions[param].values
            self.ax.arrow(0, 0, direction[0] * arrow_length, direction[1] * arrow_length,
                          head_width=0.1, head_length=0.1, fc=color, ec=color)
            self.ax.text(direction[0] * arrow_length * 1.3, direction[1]
                         * arrow_length * 1.3, param, color=color, fontsize=15)

    def plot_density(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        for strategy in self.strategy_colors:
            strategy_data = self.lda_df[self.lda_df['best_strategy'] == strategy]
            color = self.strategy_colors[strategy]
            sns.kdeplot(data=strategy_data, x='LDA Axis 1', y='LDA Axis 2', fill=True,
                        thresh=0.1, bw_adjust=1.5, alpha=0.5, ax=ax)

        # Create legend indicating strategy colors
        handles = [plt.Line2D([], [], marker='o', color='w', markerfacecolor=color, markersize=8)
                   for strategy, color in self.strategy_colors.items()]
        labels = self.strategy_colors.keys()
        ax.legend(handles, labels, loc='best')
        ax.set_xlabel('LDA Axis 1')
        ax.set_ylabel('LDA Axis 2')
        ax.set_title('Density Plot of Strategy Labels')
        self.ax = ax
        self.plot_original_axes()
