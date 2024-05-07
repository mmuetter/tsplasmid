import pandas as pd
import string
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.formula.api import ols
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt


class TwoWaySignificance:
    def __init__(self, df, group_col1, group_col2, variable, alpha=0.05):
        self.df = df
        self.group_col1 = group_col1
        self.group_col2 = group_col2
        self.variable = variable
        self.alpha = alpha
        self.anova_result = self.anova()
        if self.anova_pvalue < alpha:
            self.tukey_result = self.tukey()
            self.apply_cld()
            mean_frequencies = self.df.groupby(['strategy', 'profile'])['frequency'].mean().reset_index()
            self.cld_df.rename(columns={'group': 'strategy'}, inplace=True)
            self.cld_df = pd.merge(mean_frequencies, self.cld_df, on=['strategy', 'profile'])
            
    def anova(self):
        formula = f'{self.variable} ~ C({self.group_col1}) + C({self.group_col2}) + C({self.group_col1}):C({self.group_col2})'
        self.model = model =  ols(formula, data=self.df).fit()
        anova_results = sm.stats.anova_lm(model, typ=2)
        self.anova_pvalue = anova_results['PR(>F)'].loc["C(strategy):C(profile)"] 
        return anova_results

    def tukey(self):
        self.df['combined_group'] = self.df[self.group_col2].astype(str) + '_' + self.df[self.group_col1].astype(str)

        # Perform Tukey's HSD test on the entire dataset with the combined group identifier
        tukey_test = pairwise_tukeyhsd(endog=self.df[self.variable], 
                                       groups=self.df['combined_group'], 
                                       alpha=self.alpha)

        # Convert the results to a DataFrame
        tukey_df = pd.DataFrame(tukey_test._results_table.data[1:], columns=tukey_test._results_table.data[0])

        # Split the combined group identifier back into separate columns
        tukey_df[['profile', 'strategy1']] = tukey_df['group1'].str.split('_', expand=True)
        tukey_df[['profile2', 'strategy2']] = tukey_df['group2'].str.split('_', expand=True)

        # Filter out comparisons within the same level of strategy2
        tukey_df = tukey_df[tukey_df['profile'] == tukey_df['profile2']]

        # Drop unnecessary columns and reset index
        tukey_df.drop(columns=['group1', 'group2', 'profile2'], inplace=True)
        tukey_df.reset_index(drop=True, inplace=True)
        cols = list(tukey_df.columns)
        self.tukey_df = tukey_df[cols[-3:] + cols[:-3]]

    def write_tukey_results(self, filepath, profiles):
        if hasattr(self, 'tukey_df'):    
            df = self.tukey_df
            df = df[df.profile.isin(profiles)]
            df.to_latex(filepath)
        else:
            print("Tukey results are not available.")

    def make_cld_sub(self, df, profile):
        group1 = set(df.group1.tolist())
        group2 = set(df.group2.tolist())
        groupSet = group1 | group2
        groups = sorted(list(groupSet))

        letters = list(string.ascii_lowercase)[:len(groups)]
        cldgroups = letters

        cld = pd.DataFrame(list(zip(groups, letters, cldgroups)))
        cld[3] = ""

        for row in df.itertuples():
            if not row.reject:
                cld.iat[groups.index(df["group1"][row[0]]), 2] += cld.iat[groups.index(df["group2"][row[0]]), 1]
                cld.iat[groups.index(df["group2"][row[0]]), 2] += cld.iat[groups.index(df["group1"][row[0]]), 1]
            else:
                cld.iat[groups.index(df["group1"][row[0]]), 3] += cld.iat[groups.index(df["group2"][row[0]]), 1]
                cld.iat[groups.index(df["group2"][row[0]]), 3] += cld.iat[groups.index(df["group1"][row[0]]), 1]

        cld[2] = cld[2].apply(lambda x: "".join(sorted(set(x))))
        cld[3] = cld[3].apply(lambda x: "".join(sorted(set(x))))
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

        cld_dict = cld.set_index("groups")["labels"].to_dict()
        cld_dict = {group: {'profile': profile, 'group': group, 'letter': letter} for group, letter in cld_dict.items()}

        return cld_dict
        

    def apply_cld(self):
        cld_dicts = []
        for profile in self.df[self.group_col2].unique():
            df_subset = self.tukey_df[self.tukey_df['profile'] == profile].copy()
            df_subset.rename(columns={"strategy1": "group1", "strategy2": "group2"}, inplace=True)
            cld_dict = self.make_cld_sub(df_subset, profile)
            cld_dicts.extend(cld_dict.values())

        self.cld_df = pd.DataFrame.from_records(cld_dicts)

    def plot(self, order, include_profiles, colors, offsets, transparency = .2, fontsize_small=20, base_offset=0.03, jitter = .3, figsize=(12, 6), ax = None ):
        if ax is None:
            plt.figure(figsize=figsize)
            ax = plt.gca()

        for profile, color, offset in zip(include_profiles, colors, offsets):
            df_filtered = self.df[self.df['profile'] == profile]
            sns.boxplot(data=df_filtered, y="frequency", x="strategy", 
                        color=color, 
                        showmeans=True, meanline=True,
                        meanprops={'ls': '-', 'lw': 3, 'color': color},
                        medianprops={'visible': False},
                        whiskerprops={'visible': False},
                        showfliers=False, showbox=False, showcaps=False, 
                    order = order)

            sns.stripplot(x="strategy",
                            y="frequency",
                            data=df_filtered,
                            color=color, 
                            s=15,
                            alpha=transparency,
                            jitter = jitter,
                        order = order
                        )
            if self.anova_pvalue < self.alpha:
                self.annotate_tukey(ax, profile, color, base_offset, fontsize_small, offset)
            
            ax.set_xlabel(None)
            ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize_small)
            ax.tick_params(axis='x', labelsize=fontsize_small)
            ax.tick_params(axis='y', labelsize=fontsize_small)
            ax.set_ylim([-.1, 1.1])
        return ax

    def annotate_tukey(self, ax, profile, color, base_offset, fontsize_small, x_offset):
            profile_cld = self.cld_df[self.cld_df['profile'] == profile]
            unique_strategies = self.df['strategy'].unique()
            x_positions = dict(zip(unique_strategies, range(len(unique_strategies))))
            for _, row in profile_cld.iterrows():
                strategy = row['strategy']
                letter = row['letter']
                frequency = self.df[(self.df['profile'] == profile) & (self.df['strategy'] == strategy)]['frequency'].mean()
                y_position = frequency + base_offset
                x_position = x_positions[strategy] + x_offset
                ax.text(x_position, y_position, letter, ha='center', va='bottom', color=color, fontsize=fontsize_small)

            
    def write_anova_results(self, filepath):
        results = self.anova_result.copy()
        results.fillna("", inplace=True)
        results["PR(>F)"] = results["PR(>F)"].apply(lambda x: f"{x:.2e}" if isinstance(x, float) else x)
        results.rename(columns={
            'sum_sq': 'sum\\_sq',
            'df': 'df',
            'F': 'F',
            'PR(>F)': 'PR($>$F)'
        }, inplace=True)
        results.to_latex(filepath, index=False)
        
