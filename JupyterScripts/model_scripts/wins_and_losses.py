import matplotlib.pyplot as plt
import os
import matplotlib.patches as patches

scenario_dict = {"20220127": "preex.", "20220412": "emerg.", "20210417": "base"}


def add_bar(
    index_value, count_value, color, ax=plt, legend_label="p>0.05", bar_width=1
):
    ax.bar(index_value, count_value, color=color, label=legend_label, width=bar_width)


class WinsAndLosses:
    def __init__(
        self,
        data,
        path,
        colors,
        fontsize=16,
        axes_label_fontsize=16,
        figsize=(4, 8),
        bar_width=0.75,
        y_label_x=0,
    ):
        self.data = data
        self.colors = colors
        self.positions = data.index
        self.bar_width = bar_width
        self.axes_label_fontsize = axes_label_fontsize
        self.fontsize = fontsize
        self.path = path
        self.exp_results = False
        self.figsize = figsize
        self.strategies = list(self.data.index)
        self.y_label_x = y_label_x
        self.y_max_loss = 85
        self.y_max_win = 105

    def plot(
        self,
        title=None,
        ylabel_wins="Scenarios Won [%]",
        ylabel_losses="Scenarios Lost [%]",
        yticks=True,
        legend=True,
        ncol=2,
        hspace=0.1,
    ):
        self.fig, self.axs = plt.subplots(
            2,
            1,
            figsize=self.figsize,
            sharex=True,
            gridspec_kw={"height_ratios": [self.y_max_win, self.y_max_loss]},
        )
        self.barplot_wins(self.axs[0], ylabel_wins, title, yticks)
        self.barplot_losses(self.axs[1], ylabel_losses, yticks)
        self.fig.subplots_adjust(hspace=hspace)

        for ax in self.axs:
            self.plot_params(ax)

        if legend:
            # Creating a single legend for the entire figure
            labels = ["Single Winner", "Wins", "Single Loser", "Losses"]
            colors = [
                self.colors["single_win"],
                self.colors["total_win"],
                self.colors["single_loss"],
                self.colors["total_loss"],
            ]
            patches = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]
            self.axs[0].legend(
                patches,
                labels,
                loc="right",
                ncol=ncol,
                fontsize=self.fontsize,
                frameon=True,
            )

        # Adjusting the layout to make space for the legend if necessary

        #  plt.tight_layout(rect=[0, 0.08, 1, 0.95] if legend else [0, 0, 1, 1])
        return self.fig, self.axs

    def barplot_wins(self, ax, ylabel, title, yticks):
        df = self.data[["single_winner", "winner"]].copy()
        df.plot(
            kind="bar",
            ax=ax,
            width=self.bar_width,
            color=[self.colors["single_win"], self.colors["total_win"]],
            edgecolor="black",
            linewidth=0.5,
        )
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=self.axes_label_fontsize)
        else:
            ax.set_ylabel("", fontsize=self.axes_label_fontsize)
        ax.set_ylim(0, self.y_max_win)  # Set the upper limit for wins
        _, y_label_y = ax.yaxis.get_label().get_position()
        ax.yaxis.set_label_coords(self.y_label_x, y_label_y)
        ax.get_legend().remove()
        ax.set_title(title, fontsize=self.fontsize, fontweight="bold")
        if self.exp_results:
            self.annotate_exp(ax, "winner", self.annotate_y)

        if not yticks:
            ax.tick_params(axis="y", which="both", labelleft=False, length=0)

    def barplot_losses(self, ax, ylabel, yticks):
        df = self.data[["single_loser", "loser"]].copy()
        df = -df  # Inverting the values for downward bars
        df.plot(
            kind="bar",
            ax=ax,
            width=self.bar_width,
            color=[self.colors["single_loss"], self.colors["total_loss"]],
            edgecolor="black",
            linewidth=0.5,
        )
        ax.set_ylabel(ylabel, fontsize=self.axes_label_fontsize)
        ax.set_ylim(-self.y_max_loss, 0)  # Set the lower limit for losses
        _, y_label_y = ax.yaxis.get_label().get_position()
        ax.yaxis.set_label_coords(self.y_label_x, y_label_y)
        ax.get_legend().remove()

        if self.exp_results:
            self.annotate_exp(ax, "loser", -self.annotate_y)

        if not yticks:
            ax.tick_params(axis="y", which="both", labelleft=False, length=0)

    def annotate_exp(self, ax, case, y):
        self.get_exp_annotations(case)
        self.annotate(ax, y)

    def get_exp_annotations(self, case):
        # Collect winning scenarios for each strategy
        annotations_single = {}
        annotations_shared = {}
        for result in self.exp_results:
            strategies = result[case]
            scenario_label = scenario_dict.get(result["date"], "")
            if strategies:
                if len(strategies) > 1:
                    scenario_label += self.shared_symbol
                    for strategy in strategies:
                        annotations_shared.update({strategy: scenario_label})
                else:
                    annotations_single.update({strategies[0]: scenario_label})
        self.annotations_single = annotations_single
        self.annotations_shared = annotations_shared

    def annotate(self, ax, y, dx=0.13):
        for idx, strategy in enumerate(self.strategies):
            y_shared = y
            annotations_single = self.annotations_single.get(strategy, "")
            if ax:
                ax.annotate(
                    annotations_single,
                    xy=(idx, 0),
                    xytext=(idx - dx, y),
                    fontsize=self.fontsize * 0.9,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    color="black",
                )
            annotations_shared = self.annotations_shared.get(strategy, "")
            if ax:
                ax.annotate(
                    annotations_shared,
                    xy=(idx, 0),
                    xytext=(idx + dx, y_shared),
                    fontsize=self.fontsize * 0.9,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    color="black",
                )

    def plot_params(self, ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", which="both", length=0)
        ax.axhline(y=0, color="grey", linewidth=2, alpha=0.5, xmin=0.01, xmax=0.99)

        # Adjusting arrow coordinates and direction
        if ax.get_ylim()[0] < 0:  # For losses plot
            ax.add_patch(
                patches.FancyArrowPatch(
                    (-0.6, -self.y_max_loss),
                    (-0.6, -1),
                    arrowstyle="<-",
                    color="black",
                    mutation_scale=15,
                )
            )
        else:  # For wins plot
            ax.add_patch(
                patches.FancyArrowPatch(
                    (-0.6, 1),
                    (-0.6, self.y_max_win),
                    arrowstyle="->",
                    color="black",
                    mutation_scale=15,
                )
            )

        ax.tick_params(axis="y", labelsize=self.fontsize, rotation=0)
        ax.tick_params(axis="x", labelsize=self.fontsize, rotation=0)
        ax.set_xlabel("")

        ax.set_yticklabels([abs(int(label)) for label in ax.get_yticks()])

        # Adding horizontal grid lines
        ax.yaxis.grid(True, linestyle="-", linewidth=1, color="lightgrey", alpha=0.7)

    def add_exp_results(self, results, shared_symbol="$^+$", y=10):
        self.annotate_y = y
        self.exp_results = results
        self.shared_symbol = shared_symbol

    def save(self, path, name_add=""):
        plt.savefig(
            os.path.join(path, "wins_and_losses" + name_add + ".pdf"),
            dpi=300,
            bbox_inches="tight",  # Adjust bounding box to remove white space
            pad_inches=0.1,  # Optional: Controls padding around the figure
        )


def roman_to_int(roman):
    roman_numerals = {"I": 1, "V": 5, "X": 10, "L": 50, "C": 100, "D": 500, "M": 1000}
    total = 0
    prev_value = 0
    for numeral in reversed(roman):
        value = roman_numerals[numeral]
        if value >= prev_value:
            total += value
        else:
            total -= value
        prev_value = value
    return total


def sort_roman_numerals(roman_numerals, shared_symbol):
    def key_func(numeral):
        return roman_to_int(numeral.replace(shared_symbol, ""))

    return sorted(roman_numerals, key=key_func)


def set_title(fig, annotation, title, title_y=0.9, fontsize=20, family="serif", x=0.06):
    # Add the annotation (e.g., 'A)') at the top left of the figure
    fig.text(
        x,
        title_y,
        annotation,
        ha="left",
        va="top",
        fontsize=fontsize,
        family=family,
        weight="bold",
    )  # Adjust font size to the figure size

    # Add the title in the top center, vertically aligned with the annotation, in italic
    fig.suptitle(
        title,
        x=0.5,
        y=title_y,
        ha="center",
        va="top",
        fontsize=fontsize,
        family=family,
        fontstyle="italic",
    )  # Adjust the y position


def set_common_ylabel(fig, axs, common_label, fontsize=16, family="serif", x=0.06):
    # Remove individual y-axis labels
    for ax in axs:
        ax.set_ylabel("")

    # Add common y-axis label
    fig.text(
        x,
        0.5,
        common_label,
        va="center",
        rotation="vertical",
        fontsize=fontsize,
        family=family,
    )
