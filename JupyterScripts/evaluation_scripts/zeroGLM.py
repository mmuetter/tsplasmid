import pymc as pm
import pandas as pd
    
class BayesianZeroGLM:
    def __init__(self, df, formula, inflation_formula=None):
        self.df = df.copy()
        self.formula = formula
        self.inflation_formula = inflation_formula if inflation_formula else formula
        self.model = self.build_model()

    def build_model(self):
        with pm.Model() as model:
            # Design matrix for the count model
            X = pm.Data("X", pd.get_dummies(self.df[self.formula.split('~')[1]]))

            # Design matrix for the zero-inflation model
            X_infl = pm.Data("X_infl", pd.get_dummies(self.df[self.inflation_formula.split('~')[1]]))

            # Coefficients for the count model
            beta = pm.Normal("beta", mu=0, sigma=10, shape=X.shape[1])

            # Coefficients for the zero-inflation model
            psi = pm.Normal("psi", mu=0, sigma=10, shape=X_infl.shape[1])

            # Calculate lambda and theta
            lambda_ = pm.math.exp(pm.math.dot(X, beta))
            theta = pm.Deterministic("theta", pm.math.sigmoid(pm.math.dot(X_infl, psi)))

            # Zero-Inflated Poisson likelihood
            y_obs = pm.ZeroInflatedPoisson("y_obs", psi=theta, mu=lambda_, observed=self.df[self.formula.split('~')[0]])

        return model

    def fit(self, draws=1000, tune=1000):
        with self.model:
            self.trace = pm.sample(draws, tune=tune, return_inferencedata=True)

    def plot_trace(self):
        pm.plot_trace(self.trace)

# Example usage
# formula = 'numerator ~ pheno_group:strategy'
# inflation_formula = 'numerator ~ pheno_group' # Adjust as needed
# bayesian_glm = BayesianZeroGLM(df, formula, inflation_formula)
# bayesian_glm.fit()
# bayesian_glm.plot_trace()
