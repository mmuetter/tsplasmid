import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import string
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols

class OneWaySignificance:
    def __init__(self, df, group_col, variable, alpha=0.05):
        self.df = df
        self.group_col = group_col
        self.variable = variable
        self.alpha = alpha
        
        self.anova_result = self.anova()
        if self.anova_pvalue < alpha:
            self.tukey_result = self.tukey()
            self.tukey_df = pd.DataFrame(self.tukey_result.summary().data[1:], columns=self.tukey_result.summary().data[0])
            self.make_cld()
            
    def anova(self):
        groups = self.df[self.group_col].unique()
        V = []
        for group in groups:
            v = self.df.loc[self.df[self.group_col] == group, self.variable].values
            V.append(v)
        self.anova_return = stats.f_oneway(*V)
        f_statistic, p_value = self.anova_return


        # Calculate degrees of freedom
        df_total = len(self.df) - 1
        df_between = len(groups) - 1
        df_within = df_total - df_between

        # Calculate sum of squares
        mean_total = self.df[self.variable].mean()
        SS_total = sum((x - mean_total) ** 2 for x in self.df[self.variable])
        SS_between = sum(len(self.df[self.df[self.group_col] == group]) * (self.df[self.df[self.group_col] == group][self.variable].mean() - mean_total) ** 2 for group in groups)
        SS_within = SS_total - SS_between

        # Calculate mean squares
        MS_between = SS_between / df_between
        MS_within = SS_within / df_within

        # Create a DataFrame
        df_anova = pd.DataFrame({
            "": ["Between Groups", "Within Groups", "Total"],
            "Sum of Squares": [SS_between, SS_within, SS_total],
            "df": [df_between, df_within, df_total],
            "Mean Square": [MS_between, MS_within, ""],
            "F": [f_statistic, "", ""],
            "Sig.": [p_value, "", ""]
        })

        # Set the first column as the index without a name
        df_anova.set_index("", inplace=True)
        self.anova_pvalue = p_value
        return df_anova


    def tukey(self):
        return pairwise_tukeyhsd(endog=self.df[self.variable], groups=self.df[self.group_col], alpha=self.alpha)   


    def make_cld(self):
        df = self.tukey_df
        df["p-adj"] = df["p-adj"].astype(float)
        # Creating a list of the different treatment groups from Tukey's
        group1 = set(df.group1.tolist())  # Dropping duplicates by creating a set
        group2 = set(df.group2.tolist())  # Dropping duplicates by creating a set
        groupSet = group1 | group2  # Set operation that creates a union of 2 sets
        groups = sorted(list(groupSet))

        # Creating lists of letters that will be assigned to treatment groups
        letters = list(string.ascii_lowercase)[:len(groups)]
        cldgroups = letters

        # the following algoritm is a simplification of the classical cld,

        cld = pd.DataFrame(list(zip(groups, letters, cldgroups)))
        cld[3] = ""

        for row in df.itertuples():
            if df["p-adj"][row[0]] > (self.alpha):
                cld.iat[groups.index(df["group1"][row[0]]),
                        2] += cld.iat[groups.index(df["group2"][row[0]]), 1]
                cld.iat[groups.index(df["group2"][row[0]]),
                        2] += cld.iat[groups.index(df["group1"][row[0]]), 1]

            if df["p-adj"][row[0]] < (self.alpha):
                cld.iat[groups.index(df["group1"][row[0]]),
                        3] += cld.iat[groups.index(df["group2"][row[0]]), 1]
                cld.iat[groups.index(df["group2"][row[0]]),
                        3] += cld.iat[groups.index(df["group1"][row[0]]), 1]

        cld[2] = cld[2].apply(lambda x: "".join(sorted(x)))
        cld[3] = cld[3].apply(lambda x: "".join(sorted(x)))
        cld.rename(columns={0: "groups"}, inplace=True)

        # this part will reassign the final name to the group
        # for sure there are more elegant ways of doing this
        cld = cld.sort_values(cld.columns[2], key=lambda x: x.str.len())
        cld["labels"] = ""
        letters = list(string.ascii_lowercase)
        unique = []
        for item in cld[2]:

            for fitem in cld["labels"].unique():
                for c in range(0, len(fitem)):
                    if not set(unique).issuperset(set(fitem[c])):
                        unique.append(fitem[c])
            g = len(unique)

            for kitem in cld[1]:
                if kitem in item:
                    if cld["labels"].loc[cld[1] == kitem].iloc[0] == "":
                        cld["labels"].loc[cld[1] == kitem] += letters[g]

                    # Checking if there are forbidden pairing (proposition of solution to the imperfect script)
                    if kitem in ' '.join(cld[3][cld["labels"] == letters[g]]):
                        g = len(unique)+1

                    # Checking if columns 1 & 2 of cld share at least 1 letter
                    if len(set(cld["labels"].loc[cld[1] == kitem].iloc[0]).intersection(cld.loc[cld[2] == item, "labels"].iloc[0])) <= 0:
                        if letters[g] not in list(cld["labels"].loc[cld[1] == kitem].iloc[0]):
                            cld["labels"].loc[cld[1] == kitem] += letters[g]
                        if letters[g] not in list(cld["labels"].loc[cld[2] == item].iloc[0]):
                            cld["labels"].loc[cld[2] == item] += letters[g]

        cld = cld.sort_values("labels")
        cld.drop(columns=[1, 2, 3], inplace=True)
        self.letters = (cld.set_index("groups"))

    def plot(self, ax=None, ylim=None, transparancy=0.2, fontsize_small=20, figsize=(12, 8), bar_color='red', dot_color='blue', base_offset=0.05, order_col="strategy", order=False, jitter = .2, violin = True, zorder = 1, strip = False, dx = 0):
        if ax is None:
            plt.figure(figsize=figsize)
            ax = plt.gca()

        if order:
            # Creating a copy of the DataFrame to avoid SettingWithCopyWarning
            df = self.df.copy()
            df[order_col] = pd.Categorical(df[order_col], categories=order, ordered=True)
            df.sort_values(order_col, inplace=True)
        else:
            df = self.df

        means = self.df[[self.group_col, self.variable]].groupby([self.group_col]).mean()

        if violin:
            sns.violinplot(x=self.group_col, 
                y=self.variable, 
                data=df,
                inner=None, 
                ax=ax)

        for violin in ax.collections:
            violin.set_edgecolor(dot_color)
            violin.set_alpha((transparancy+1)/2)
            violin.set_facecolor('none')

        sns.boxplot(showmeans=True,
                    meanline=True,
                    meanprops={'color': bar_color, 'ls': '-', 'lw': 3},
                    medianprops={'visible': False},
                    whiskerprops={'visible': False},
                    x=self.group_col,
                    y=self.variable,
                    data=df,
                    showfliers=False,
                    showbox=False,
                    showcaps=False,
                    ax=ax,
                    color=bar_color)

        if strip:
            sns.stripplot(x=self.group_col,
                        y=self.variable,
                        data=df,
                        color=dot_color,
                        s=15,
                        alpha=transparancy,
                        ax=ax, 
                        jitter = jitter,
                        zorder = zorder)

        if ylim:
            ax.set_ylim(ylim)

        if self.anova_pvalue < self.alpha:
            self.add_tukey(ax, means, bar_color, base_offset*df[self.variable].max(), fontsize_small, x_offset = dx)
        return ax
    


    def add_tukey(self, ax, means, color, base_offset,fontsize, x_offset):
        if hasattr(self, 'letters'):
            for xtick in ax.get_xticks():
                label = ax.get_xticklabels()[xtick].get_text()
                pos_vec = means[self.variable] - 0.5
                offset = base_offset if pos_vec[label] else -base_offset
                string = self.letters.loc[label, "labels"]
                ax.text(xtick + x_offset, means.loc[label, self.variable] + offset, string,
                        horizontalalignment='center', color=color, weight='semibold', fontsize=fontsize)

    def write_anova_results(self, filepath, decimal_places=3):
        modified_df = self.anova_result.applymap(lambda x: '$< 0.001$' if isinstance(x, float) and x < 0.001 else x)
        float_format_func = f"{{:0.{decimal_places}f}}".format
        modified_df.to_latex(filepath, index=True, float_format=float_format_func)

    def write_tukey_results(self, filepath, decimal_places=3):
        self.tukey_df.to_latex(filepath, index=False, float_format=f"{{:0.{decimal_places}f}}".format)
