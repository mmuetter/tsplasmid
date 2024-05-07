import pandas as pd



class Superconjugation:
    def __init__(self, experiment, fontsize = 21, tstart = 4):
        self.fontsize = fontsize
        self.experiment = experiment
        data = experiment.data.rename(
            columns={"transfered_phenotype": "transferred_phenotype"})
        data = data[data.transfer_n >= tstart]
        data["encounter"] = data.apply(lambda x:
                                       ("A_r" in x.added_strains) &
                                       ("B_r" in x.added_strains) &
                                       ("AB_r" not in x.added_strains) &
                                       ("A&B" not in x.added_strains) &
                                       ("Other" not in x.added_strains), axis=1)
        encounter = data[data.encounter].copy()
        encounter["emergence"] = encounter.phenotype == "AB_r"
        self.encounter = encounter

    # select wells were A-r and B-r encounter
    # exclude wells with preexisting double resistance: 𝐴𝐵𝑟
    # also exclude wells were it is uncertain if there is a preexisting double resistance
    # Find the probability of encounter between Ar and Br
    # Devide by 376 patients per transfer

    def encounter_frequencies(self, alpha = 0.05, lim = (0, .15)):
        cols = ["strategy", "transfer_n", "encounter"]
        summary = self.encounter[cols].groupby(cols[:-1]).count()
        summary = summary.reset_index()
        summary["first"] = summary.transfer_n == 1
        summary["f"] = summary.encounter/376
        self.encounter_summary = summary
        self.superinfections = summary[summary["first"] == False]
        self.superinfections_first = summary[summary["first"]]
  
   
    def emergence_per_encounter(self):
        cols = ["treatment_with", "transfer_n", "strategy",  "emergence"]
        summary_ab = self.encounter[cols].groupby(cols[:-1]).sum()
        summary_ab["encounters"] = self.encounter[cols].groupby(
            cols[:-1]).count()
        summary_ab = summary_ab.reset_index()
        summary_ab["f"] = summary_ab.emergence / summary_ab.encounters
        self.emergence_summary = summary_ab


    def calculate_expected_frequencies(self):
        cols = ["strategy", "f"]
        enc = self.encounter_summary[cols].groupby(cols[:-1]).mean()
        cols = ["treatment_with", "f"]
        em = self.emergence_summary[cols].groupby(cols[:-1]).mean()

        treatment = []
        data = self.experiment.data
        cyc = data[data.strategy == "Cycling"]
        for t in cyc.transfer_n.unique():
            treatment.append(
                cyc[cyc.transfer_n == t].treatment_with.unique()[0])
        n = len(treatment)
        n_A = treatment.count("A")
        n_B = treatment.count("B")
        treatment

        expectations = {
            "No treatment": enc.loc["No treatment", "f"]*em.loc["none", "f"],
            "Mono A": enc.loc["Mono A", "f"]*em.loc["A", "f"],
            "Mono B": enc.loc["Mono B", "f"]*em.loc["B", "f"],
            "Mixing": enc.loc["Mixing", "f"]*(em.loc["A", "f"] + em.loc["B", "f"])/2,
            "Combination": enc.loc["Combination", "f"]*em.loc["AB", "f"],
            "Cycling": enc.loc["Cycling", "f"]*(n_A * em.loc["A", "f"] + n_B * em.loc["B", "f"])/n
        }
        df = pd.DataFrame.from_dict(
            expectations, columns=["f"], orient="index")
        self.expectations = df
