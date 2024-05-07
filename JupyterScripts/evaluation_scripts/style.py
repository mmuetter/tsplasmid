
from matplotlib import pyplot as plt
import seaborn as sns


class Style:
    def __init__(self, theme="white", fontsize=12, figsize=(12, 8), palette="pastel"):
        self.theme = theme
        self.fontsize = fontsize
        self.figsize = figsize
        self.palette = palette
        self.update_rcparams()
        self.create_figure_with_style(dummy = True)
    
    def update_rcparams(self):
        self.set_theme(self.theme)
        self.set_fontsize(self.fontsize)
        self.set_figsize(self.figsize)
        self.set_palette(self.palette)

    def set_theme(self, theme):
        if theme not in ["white", "black"]:
            raise ValueError("Invalid theme. Use 'white' or 'black'.")
        self.theme = theme
        self.set_color_params()

    def set_fontsize(self, fontsize):
        if not isinstance(fontsize, int) or fontsize <= 0:
            raise ValueError("Fontsize should be a positive integer.")
        self.fontsize = fontsize
        plt.rcParams["font.size"] = self.fontsize

    def set_figsize(self, figsize):
        if (
            not isinstance(figsize, tuple)
            or len(figsize) != 2
            or not all(isinstance(val, (int, float)) and val > 0 for val in figsize)
        ):
            raise ValueError("Figsize should be a tuple of two positive numbers.")
        self.figsize = figsize
        plt.rcParams["figure.figsize"] = self.figsize

    def set_palette(self, palette):
        self.palette = palette
        if self.palette == "pastel":
            sns.set_palette("pastel")
        else:
            sns.set_palette(palette)

    def set_color_params(self):
        if self.theme == "white":
            plt.rcParams.update(
                {
                    "text.color": "black",
                    "axes.labelcolor": "black",
                    "axes.edgecolor": "black",
                    "axes.facecolor": "white",
                    "axes.titlecolor": "black",
                    "xtick.color": "black",
                    "ytick.color": "black",
                }
            )
        elif self.theme == "black":
            plt.rcParams.update(
                {
                    "text.color": "white",
                    "axes.labelcolor": "white",
                    "axes.edgecolor": "white",
                    "axes.facecolor": "black",
                    "axes.titlecolor": "white",
                    "xtick.color": "white",
                    "ytick.color": "white",
                }
            )

    def create_figure_with_style(
        self, title=None, xlabel=None, ylabel=None, figsize=False, dummy = False
    ):
        self.update_rcparams()
        if not figsize:
            figsize = self.figsize

        fig, ax = plt.subplots(figsize=figsize)
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if dummy:
            plt.clf()
        else: 
            return fig, ax
    
