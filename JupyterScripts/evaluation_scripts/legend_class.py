import pandas as pd
import matplotlib.lines as mlines

class Legend:
    def __init__(self, fig, ax, aesthetics):
        self.fig = fig
        self.ax = ax
        self.aesthetics = aesthetics
        self.handles, self.labels = self.ax.get_legend_handles_labels()
        self.label_dict = {
            "U":"$U$",
            "S":"$S$",
            "A_r":"$A_r$",
            "B_r":"$B_r$",
            "A&B":"$(A_r&B_r)$",
            "AB_r":"$AB_r$",
            "Other":r'$\mathit{Other}$'
        }

    def sort_handles_labels(self):
        rank = list(range(len(self.aesthetics.legend_order)))
        order_dict = dict(zip(self.aesthetics.legend_order, rank))
        tmp = pd.DataFrame()
        tmp["handles"] = self.handles
        tmp["names"] = self.labels
        tmp["rank"] = tmp.apply(lambda x: order_dict[x["names"]], axis=1)
        tmp.sort_values("rank", inplace=True)
        tmp["labels"] = tmp.names.apply(lambda x: self.make_labels(x))
        self.handles = list(tmp.handles)
        self.labels = list(tmp.labels)

    def make_labels(self, x):
        if x in self.label_dict.keys():
            return self.label_dict[x]
        else:
            return x

    def show(self):
        self.sort_handles_labels()
        legend_font_props = {'size': self.aesthetics.legend_size}
        self.fig.legend(self.handles, self.labels, 
                        markerscale=self.aesthetics.legend_scale, 
                        bbox_to_anchor=(0.5, -0.1), 
                        loc='lower center', 
                        ncol=len(self.labels), 
                        prop=legend_font_props, 
                        frameon=False) 
    #  Plot Phyno
    def make_exp_data_entry(self, s=2):
        info = self.aesthetics.legend["exp"]
        entry = info["entry"]
        s = info["s"]
        color = info["color"]
        p1, = self.ax.plot([], [], c=color, marker="|", markersize=5*s,
                           linestyle='None')  # notice the comma!
        p2, = self.ax.plot([], [], c=color, marker="o", markersize=s,
                           linestyle='None')  # notice the comma!
        handle = ((p1, p2),)
        label = entry
        self.handles.append(handle)
        self.labels.append(label)

    def make_val_legend(self):
        info = self.aesthetics.legend["val"]
        entry = info["entry"]
        linewidth = info["linewidth"]
        color = info["color"]
        h, = self.ax.plot([], [], color=color,
                          linestyle="dashdot", label=entry, linewidth=linewidth)
        self.handles.append(h)
        self.labels.append(entry)

    def make_var_legend(self):
        info = self.aesthetics.legend["var"]
        entry = info["entry"]
        linewidth = info["linewidth"]  # Controls thickness of p1
        color = info["color"]
        blockwidth = info["blockwidth"]  # Controls thickness of p2

        if color == "black":
            alpha = 0.1
        else:
            alpha = 0.4

        # p1 - Adjusting linewidth or linestyle could change perceived length
        p1_line = mlines.Line2D([], [], c='black', linestyle="-", linewidth=linewidth)

        # p2 - Adjusting linewidth changes thickness, alpha for transparency
        p2_line = mlines.Line2D([], [], c=color, linestyle=(0, (5, 10)), linewidth=blockwidth, alpha=alpha)

        custom_handle = (p1_line, p2_line)

        self.handles.append(custom_handle)
        self.labels.append(entry)