import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson, gamma, binom, norm
from itertools import product
from evaluation_scripts.cld import resolve_significance_table
import itertools
import os

class GLM:
    def __init__(self, experiment, pheno_group, dt = 4):
        self.pheno_group = pheno_group
        self.experiment = experiment
        self.dt = dt
        self.summarize_data()  
              


    def summarize_data(self):
        pheno_dict = {"S":"other", 
            "Other":"other",
            "A_r":"single resistant", 
            "B_r":"single resistant", 
            "A&B":"single resistant", 
            "U":"uninfected", 
            "AB_r": "double resistant"}
        exp = self.experiment
        Data = exp.data
        Data["pheno_group"] = Data.phenotype.map(pheno_dict)
        Data = exp.data[exp.data.pheno_group == self.pheno_group]
        t_max = max(Data.transfer_n)
        data = Data[Data.transfer_n > t_max-self.dt]
        collect = []
        for s, r in product(exp.strategies, [0, 1, 2, 3]):
            l = len(data[(data.strategy == s) & (data.rep == r)])
            collect.append({
                "strategy":s,
                "rep":r,
                "count":l,
                "count_mean":l/self.dt,
                "patients":self.dt*94,
                "f" : l/(self.dt*94)
            })
        self.df = pd.DataFrame().from_records(collect)

    
    def fit_glm(self, formula, family, link=None, fit_disp=False, intercept:tuple = ("strategy", "No treatment"), offset = 0):
        self.sort(intercept)
        self.formula = formula
        self.predictors = self.formula.split('~')[1]
        self.variable = self.formula.split('~')[0]
        self.mean_variance = self.calculate_mean_variance()
        #link_instance = link() if link is not None else None

        try:
            self.family  = getattr(sm.families, family)
            self.family = self.family(link=link)
        except AttributeError:
            raise ValueError(f"Unsupported family: {family}")

        if fit_disp:
            self.theta = self.optimize_theta()
            self.family.variance = lambda mu: self.theta * mu
        
        self.df_offset = self.df.copy()
        self.df_offset[self.variable] = self.df[self.variable] + offset
        self.model = smf.glm(formula=self.formula, data=self.df_offset, family=self.family).fit()
        self.null_model = self.get_null_model() 
        self.make_summary()       
        self.format_summaries()
        self.check_significance()
        self.transform_significance()
        self.cld(self.summaries[1])
        self.add_means()

    def make_summary(self):
        self.summaries = []
        for summary in self.model.summary().tables:
            self.summaries.append(pd.DataFrame(summary.data[1:], columns=summary.data[0]))
        
        summary = self.summaries[1]
        summary.loc[:, ["coef", "std err", "[0.025", "0.975]"]] = summary.loc[:, ["coef", "std err", "[0.025", "0.975]"]].astype(float)

        intercept = summary.iloc[0]
        ref = {
            "":"strategy[T.No treatment] (reference)",
            "coef" : 0, 
            "std err":intercept["std err"], 
            "z":"", 
            "P>|z|":"", 
            "[0.025":intercept["[0.025"]-intercept.coef,
            "0.975]":intercept["0.975]"]-intercept.coef }
        summary.loc[6, :] = ref
        self.summaries[1] = add_strategy_column(summary)
    
    def add_means(self):
        df = self.df.groupby(self.predictors).mean()
        summary =  self.summaries[1]
        summary[self.variable + "_mean"] = summary[self.predictors].map(df[self.variable]/self.dt)
        self.summaries[1] =  summary

    def get_null_model(self):
        null_formula = f'{self.variable} ~ 1'
        null_model = smf.glm(formula=null_formula, data=self.df_offset, family=self.family).fit()
        return null_model

    def sort(self, intercept:tuple):
        intercept = ("strategy", "No treatment")
        entries = set(self.df[intercept[0]])
        entries.remove(intercept[1])
        self.df['strategy'] = pd.Categorical(self.df['strategy'], categories=[intercept[1]] + list(entries), ordered=True)


    def get_variance(self, result, var):
        subgroup_means = result.groupby(["exp", "strategy", "pheno_group"]).mean()[var].reset_index(name='subgroup_mean')
        result_with_means = pd.merge(result, subgroup_means, on=["exp", "strategy", "pheno_group"])
        result_with_means['residuals'] = result_with_means[var] - result_with_means['subgroup_mean']
        subgroup_variances = result.groupby(["exp", "strategy", "pheno_group"])[var].var().reset_index(name='subgroup_variance')
        result_with_variances = pd.merge(result_with_means, subgroup_variances, on=["exp", "strategy", "pheno_group"])
        return result_with_variances.dropna()       

    def plot_histograms(self):
        # Create a copy of the dataframe
        df_copy = self.df.copy()

        # Splitting the predictors string and stripping whitespace
        predictors = [pred.strip() for pred in self.predictors.split(':')]

        # Create a new column 'subgroup' by concatenating the predictor columns
        df_copy['subgroup'] = df_copy[predictors].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

        # Plotting histograms using seaborn's FacetGrid
        g = sns.FacetGrid(df_copy, col="subgroup", col_wrap=3, sharex=False, sharey=False)
        g.map(sns.histplot, self.variable)
        plt.show()
    
    def plot_distribution(self, ax=None):
        # Extract the observed data
        data = self.df[self.variable]

        # Create a new figure and axis if ax is not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.figure
        sns.histplot(data, kde=False, color='blue', stat='density', bins=30, label='Observed Data', ax=ax)

        # Range for x-axis
        x = np.linspace(min(data), max(data), 100)

        # Plotting the fitted distribution
        if isinstance(self.family, sm.families.Poisson):
            # For Poisson, the mean parameter is the rate (lambda)
            rate = self.model.params[0]
            y = poisson.pmf(np.floor(x), mu=rate)
            ax.plot(x, y, 'r-', label='Fitted Poisson Distribution')
        elif isinstance(self.family, sm.families.Gamma):
            # For Gamma, parameters might include shape and scale
            shape, scale = self.model.params[0], 1  # Adjust as needed
            y = gamma.pdf(x, a=shape, scale=scale)
            ax.plot(x, y, 'r-', label='Fitted Gamma Distribution')
        elif isinstance(self.family, sm.families.Binomial):
            # For Binomial, the parameter is the probability of success (p)
            p = self.model.params[0]
            y = binom.pmf(np.floor(x), n=1, p=p)
            ax.plot(x, y, 'r-', label='Fitted Binomial Distribution')
        elif isinstance(self.family, sm.families.Gaussian):
            # For Gaussian, parameters are mean and standard deviation
            mean, std = self.model.params[0], np.sqrt(self.model.scale)
            y = norm.pdf(x, loc=mean, scale=std)
            ax.plot(x, y, 'r-', label='Fitted Gaussian Distribution')
        # Add other distributions here if needed

        # Labels and legend
        ax.set_xlabel(self.variable)
        ax.set_ylabel('Density')
        ax.set_title(f'Distribution of {self.variable} with Fitted Model Overlay')
        ax.legend()



    def calculate_mean_variance(self):
        predictor_names = [pred.strip() for pred in self.predictors.replace('+', ':').split(':')]
        grouped = self.df.groupby(predictor_names)
        mean_variance = grouped[self.variable].agg(['mean', 'var'])
        return mean_variance
        
    def optimize_theta(self):
        def to_minimize(theta):
            predicted_variance = theta * self.mean_variance['mean']
            return ((self.mean_variance['var'] - predicted_variance) ** 2).sum()
        result = minimize(to_minimize, x0=1.0, method='Nelder-Mead')
        return result.x[0]


    def report_theta(self):
        # Reporting and plotting for theta
        print(f"Optimized Theta: {self.theta}")

        # Plot observed vs predicted variance
        predicted_variance = self.theta * self.mean_variance['mean']

        plt.scatter(self.mean_variance['mean'], self.mean_variance['var'], label='Observed Variance')
        plt.plot(self.mean_variance['mean'], predicted_variance, color='red', label='Predicted Variance')
        plt.xlabel('Mean')
        plt.ylabel('Variance')
        plt.title('Mean vs Variance with Fitted Theta')
        plt.legend()
        plt.show()

    def qq_plot(self, ax=None, figsize = (12,8)):
        if ax is None:
            plt.figure(figsize=figsize)
            ax = plt.gca()
        residuals = self.model.resid_pearson
        sm.qqplot(residuals, line='45', ax=ax)
        ax.set_title('Q-Q Plot of Residuals')


    def residual_plot(self, figsize=(12, 6)):
        fitted_values = self.model.fittedvalues
        residuals = self.model.resid_pearson
        predictor = self.formula.split('~')[1].strip()

        plt.figure(figsize=figsize)
        if predictor in self.df.select_dtypes(include=['category', 'object']).columns:
            sns.boxplot(x=self.df[predictor], y=residuals)
            plt.axhline(y=0, color='red', linestyle='--')
            plt.title('Residuals by Category')
        else:
            plt.scatter(fitted_values, residuals)
            plt.xlabel('Fitted Values')
            plt.ylabel('Residuals')
            plt.axhline(y=0, color='red', linestyle='--')
            plt.title('Residuals vs Fitted')
        plt.show()

    def categorical_plot(self, figsize=(12, 6)):
        predictor = self.formula.split('~')[1].strip()
        if predictor in self.df.select_dtypes(include=['category', 'object']).columns:
            fitted_values = self.model.fittedvalues
            categories = self.df[predictor]

            plt.figure(figsize=figsize)
            sns.stripplot(x=categories, y=self.df[self.formula.split('~')[0]], jitter=True, alpha=0.5)
            sns.pointplot(x=categories, y=fitted_values, color='red', join=False)
            plt.ylabel('Fitted Values')
            plt.title('Actual vs. Fitted Values for Categorical Predictor')
            plt.show()

    def mean_variance_plot(self, figsize=(12, 6)):
        plt.figure(figsize=figsize)
        plt.scatter(self.mean_variance['mean'], self.mean_variance['var'], label='Observed Variance')

        if isinstance(self.family, sm.families.Gaussian):
            constant_variance = np.full_like(self.mean_variance['mean'], self.model.scale)
            plt.plot(self.mean_variance['mean'], constant_variance, color='red', label='Predicted (Constant Variance)')
        elif isinstance(self.family, sm.families.Poisson):
            plt.plot(self.mean_variance['mean'], self.mean_variance['mean'], color='red', label='Predicted (Variance = Mean)')
        elif 'theta' in dir(self):  # Check if theta is defined for Quasi-Poisson
            predicted_variance = self.theta * self.mean_variance['mean']
            plt.plot(self.mean_variance['mean'], predicted_variance, color='red', label='Predicted (Quasi-Poisson Variance)')
        elif isinstance(self.family, sm.families.Binomial):
            p = self.model.predict()  # Predicted probabilities
            individual_variance = p * (1 - p)
            predictor = self.formula.split('~')[1].strip()
            df_with_preds = self.df.copy()
            df_with_preds['predicted_variance'] = individual_variance
            predicted_variance_by_group = df_with_preds.groupby(predictor)['predicted_variance'].mean()
            plt.plot(self.mean_variance['mean'], predicted_variance_by_group, color='red', label='Predicted Variance (Binomial)')
        plt.xlabel('Mean')
        plt.ylabel('Variance')
        plt.title('Mean vs. Variance')
        plt.legend()
        plt.show()

    def plot_summary(self, figsize=(12, 6)):
        self.residual_plot(figsize)
        self.qq_plot(figsize)
        self.mean_variance_plot(figsize)
        predictor = self.formula.split('~')[1].strip()
        if predictor in self.df.select_dtypes(include=['category', 'object']).columns:
            self.categorical_plot(figsize)

    def check_significance(self):
        summary = self.summaries[1]
        summary = summary[summary[""] != "Intercept"]
        groups = summary.iloc[:, 0].tolist()
        summary = summary.set_index("")
        permutations = [(x, y) for x, y in itertools.permutations(groups, 2) if x != y]
        significance_results = []
        
        for group_a, group_b in permutations:
            mean_a = summary.loc[group_a, 'coef']
            ci_lower_b, ci_upper_b = summary.loc[group_b, '[0.025'], summary.loc[group_b, '0.975]']
            is_significant = not (ci_lower_b <= mean_a <= ci_upper_b)
            significance_results.append((group_a, group_b, mean_a, ci_lower_b, ci_upper_b, is_significant))

        significance_df = pd.DataFrame(significance_results, columns=['Group A', 'Group B', 'Mean A', 'CI Lower B', 'CI Upper B', 'Significant'])
        self.significance_df = significance_df

    def transform_significance(self):
        # Make a copy of the original significance DataFrame
        significance_df_copy = self.significance_df.copy()

        # Initialize a list to store transformed results
        transformed_records = []

        # Iterate through each row in the copied significance DataFrame
        for index, row in significance_df_copy.iterrows():
            group1, group2 = row['Group A'], row['Group B']
            sig1 = row['Significant']

            # Find the corresponding reverse pair in the copy
            reverse_pair = significance_df_copy[(significance_df_copy['Group A'] == group2) & (significance_df_copy['Group B'] == group1)]
            if not reverse_pair.empty:
                sig2 = reverse_pair['Significant'].values[0]
                # Add to the list of transformed records
                transformed_records.append({'Group 1': group1, 'Group 2': group2, 'Sig Group1 -> Group2': sig1, 'Sig Group2 -> Group1': sig2})

        # Create a DataFrame from the list of records
        self.transformed_significance = pd.DataFrame.from_records(transformed_records)
        self.transformed_significance["significance"] = self.transformed_significance["Sig Group1 -> Group2"] & self.transformed_significance["Sig Group2 -> Group1"]
        self.significance_pivot = self.transformed_significance.pivot(index = "Group 1", columns = "Group 2", values = "significance")


    def cld(self, r):
        L = ["a", "b", "c", "d", "e", "f"]
        r["cld"] = ""
        i = r[""][1:]
        M = self.significance_pivot.loc[i, i]
        groups = {frozenset(s) for s in resolve_significance_table(M)} 
        groups = [set(s) for s in groups]
        for i, g in enumerate(groups):
            idx = r[r[""].isin(g)].index
            r.loc[idx, 'cld'] = r.loc[idx, "cld"].apply(lambda x: x + L[i])
        return r

    def plot_transformed_frequencies(self, figsize=(12, 6)):
            # Apply the link function to transform the frequencies
            transformed_f = self.model.family.link.inverse(self.df[self.variable])

            # Determine the label for the y-axis based on the link function
            ylabel = f'{type(self.model.family.link).__name__}(f)'

            # Create the plot
            plt.figure(figsize=figsize)
            sns.stripplot(x=self.df['strategy'], y=transformed_f, jitter=True, alpha=0.5)

            # Overlay with boxplot only showing the mean
            sns.boxplot(x=self.df['strategy'], y=transformed_f, showcaps=False, boxprops={'facecolor':'None'}, 
                        showfliers=False, whiskerprops={'linewidth':0}, meanline=True, meanprops={"color": "red"})

            plt.ylabel(ylabel)
            plt.title('Transformed Frequencies by Strategy')

    def report(self):
        path = os.path.join(self.experiment.pathes["tables"], self.experiment.exp + "_glm.tex")
        summary = self.summaries[1].iloc[0:-1, 0:-1]
        summary.to_latex(path, float_format="%.3f", index=False)

    def format_summaries(self):
        formatted_summaries = []
        for summary in self.summaries:
            # Round float columns to three decimal places
            for col in summary.select_dtypes(include=['float']):
                summary[col] = summary[col].round(3)

            # Enclose column names containing '>' in '$'
            summary.columns = [f"${col}$" if '>' in col else col for col in summary.columns]

            formatted_summaries.append(summary)
        
        return formatted_summaries
    
    def plot_strategies(self, color = "cyan", ycol = "count_mean", xcol = "strategy", transparency = .2, fontsize_small=20, base_offset=0.03, x_offset = 0, jitter = .3, figsize=(16, 6), ax = None , ylim = (-1, 94), annotate = True):
        if ax is None:
                plt.figure(figsize=figsize)
                ax = plt.gca()
        
        df = self.df
        sns.boxplot(data=df, y=ycol, x=xcol, 
                    color=color, 
                    showmeans=True, meanline=True,
                    meanprops={'ls': '-', 'lw': 3, 'color': color},
                    medianprops={'visible': False},
                    whiskerprops={'visible': False},
                    showfliers=False, showbox=False, showcaps=False, 
                order = self.experiment.strategies)

        sns.stripplot(x=xcol,
                        y=ycol,
                        data=df,
                        color=color, 
                        s=15,
                        alpha=transparency,
                        jitter = jitter,
                    order = self.experiment.strategies
                    )
        
        ax.set_xlabel(None)
        ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize_small)
        ax.tick_params(axis='x', labelsize=fontsize_small)
        ax.tick_params(axis='y', labelsize=fontsize_small)
        ax.set_ylim(*ylim)

        if annotate:
            ax = self.annotate(ax, color, base_offset, x_offset, fontsize_small)
        return ax
    
    def annotate(self, ax, color, base_offset, x_offset, fontsize_small):
        summary = self.summaries[1].copy()
        summary[summary[self.predictors].isnull() == False]
        summary = summary.set_index(self.predictors)
        x_labels = [label.get_text() for label in ax.get_xticklabels()]        
        for i, predictor in enumerate(x_labels):
                letter = summary.loc[predictor, 'cld']
                mean = summary.loc[predictor, self.variable + '_mean']
                y_position = mean + base_offset
                x_position = i + x_offset
                ax.text(x_position, y_position, letter, ha='center', va='bottom', color=color, fontsize=fontsize_small)
        return ax

def add_strategy_column(df):
            # Extract strategy names from the index or a specific column
            df['strategy'] = df[""].str.extract(r'\[T\.(.*?)\]')[0]
            return df