import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import minimize

class multiGLM:
    def __init__(self, df, formula, family, link=None, fit_disp=False):
        df['strategy'] = pd.Categorical(df['strategy'], categories=['No treatment', 'Combination', 'Cycling', 'Mixing', 'Mono A', 'Mono B'], ordered=True)
        self.df = df
        self.formula = formula
        self.mean_variance = self.calculate_mean_variance()
        link_instance = link() if link is not None else None

        try:
            family_class = getattr(sm.families, family)
            self.family = family_class(link=link_instance)
        except AttributeError:
            raise ValueError(f"Unsupported family: {family}")

        if fit_disp:
            self.theta = self.optimize_theta()
            self.family.variance = lambda mu: self.theta * mu
           
        self.model = smf.glm(formula=self.formula, data=self.df, family=self.family).fit()
    
    def optimize_theta(self):
        def to_minimize(theta):
            predicted_variance = theta * self.mean_variance['mean']
            return ((self.mean_variance['var'] - predicted_variance) ** 2).sum()
        result = minimize(to_minimize, x0=1.0, method='Nelder-Mead')
        return result.x[0]

    def calculate_mean_variance(self):
        grouped = self.df.groupby(self.df[self.formula.split('~')[1].strip()])
        mean_variance = grouped['f'].agg(['mean', 'var'])
        return mean_variance

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

    def calc_confidence_intervals(self):
        summary_data = self.model.summary().tables[1].data
        summary_df = pd.DataFrame(summary_data[1:], columns=summary_data[0])
        summary_df.iloc[:, 1:] = summary_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
        conf = self.model.conf_int()
        conf.columns = ['CI Lower', 'CI Upper']
        conf.index = summary_df.index  
        self.summary = summary_df.join(conf)

    def qq_plot(self, figsize=(12, 6)):
        residuals = self.model.resid_pearson
        plt.figure(figsize=figsize)
        sm.qqplot(residuals, line='45')
        plt.title('Q-Q Plot of Residuals')
        plt.show()

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

    def cld(self):
        # Ensure confidence intervals are calculated
        if not hasattr(self, 'summary'):
            self.calc_confidence_intervals()

        strategies = self.summary.index[1:]  # Exclude intercept
        conf_intervals = self.summary[['CI Lower', 'CI Upper']].iloc[1:]

        strategy_letters = {strategy: '' for strategy in strategies}
        assigned_letters = set()

        for i, strategy1 in enumerate(strategies):
            if strategy_letters[strategy1] == '':
                # Assign a new letter if not already assigned
                current_letter = 'a'
                while current_letter in assigned_letters:
                    current_letter = chr(ord(current_letter) + 1)

                strategy_letters[strategy1] = current_letter
                assigned_letters.add(current_letter)

            for j, strategy2 in enumerate(strategies):
                if j > i:
                    overlap = not (conf_intervals.iloc[i]['CI Upper'] < conf_intervals.iloc[j]['CI Lower'] or
                                   conf_intervals.iloc[j]['CI Upper'] < conf_intervals.iloc[i]['CI Lower'])

                    if overlap:
                        # Share the letter if overlapping and not already assigned
                        letter_to_assign = strategy_letters[strategy1]
                        if letter_to_assign not in strategy_letters[strategy2]:
                            strategy_letters[strategy2] += letter_to_assign

        self.summary['cld'] = self.summary.index.map(strategy_letters)

    def plot_transformed_frequencies(self, figsize=(12, 6)):
        # Apply the link function to transform the frequencies
        transformed_f = self.model.family.link.inverse(self.df['f'])

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
    