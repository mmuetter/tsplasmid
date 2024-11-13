from evaluation_scripts.base import load_toml, removekey
from matplotlib import pyplot as plt
import matplotlib as mpl


class Aesthetics:
    def __init__(
        self,
        style="white",
        title=False,
        fontsize=30,
        box=False,
        grid_alpha=0.2,
        x_tick_space=2,
        n=1,
        legend_scale=3,
        legend_size=20,
        legend_fontsize=15,
        grid=True,
        font_family="sans-serif",
        sans_serif_font="Helvetica",
        math_fontset="stix",
    ):
        self.title = title
        self.style = style
        self.box = box
        self.fontsize_small = 0.8 * fontsize
        self.fontsize_large = 1.4 * fontsize
        self.fontsize = fontsize
        self.legend_fontsize = legend_fontsize
        self.grid_alpha = grid_alpha
        self.x_tick_space = x_tick_space
        self.plot_n_times = n
        self.phenotype_colors = load_toml("phenotype_colors.toml")
        self.grid = grid
        self.legend = {}
        self.phenotype_order = ["U", "S", "A_r", "B_r", "A&B", "AB_r", "Other"]
        self.legend_order = self.phenotype_order.copy()
        self.legend_scale = legend_scale
        self.legend_size = legend_size

        # Font configuration
        self.font_family = font_family
        self.sans_serif_font = sans_serif_font
        self.math_fontset = math_fontset

        self.set_rcparams()
        _ = plt.draw()
        _ = plt.clf()

    def time_series(self):
        self.axs_dict = {
            "No treatment": (0, 0),
            "Mono B": (0, 1),
            "Mixing": (0, 2),
            "Mono A": (1, 0),
            "Cycling": (1, 1),
            "Combo": (1, 2),
        }
        self.xlabel = "transfer"
        self.ylabel = "proportion of resistance profiles"

    def add_legend(self, info, case):
        info.update({"color": self.rcparams["lines.color"]})
        self.legend.update({case: info})
        self.legend_order.append(info["entry"])

    def set_rcparams(self):
        black = plt.style.library["dark_background"].copy()
        white = black.copy()
        for key in black:
            if black[key] == "black":
                white.update({key: "white"})
            elif black[key] == "white":
                white.update({key: "black"})
        if self.style == "white":
            rcparams = white
        else:
            rcparams = black

        # Apply custom font settings
        rcparams.update(
            {
                "font.size": self.fontsize,
                "font.family": self.font_family,
                "font.sans-serif": [self.sans_serif_font],
                "mathtext.fontset": self.math_fontset,
            }
        )

        if not self.box:
            rcparams["axes.edgecolor"] = rcparams["axes.facecolor"]

        rcparams = removekey(rcparams, "grid.color")
        self.rcparams = rcparams
        mpl.rcParams.update(**self.rcparams)
