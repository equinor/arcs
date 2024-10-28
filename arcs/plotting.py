import contextlib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

import numpy as np
import matplotlib.colors as colors
import os
import matplotlib


def plot_difference(
    axis,
    mean_data,
    temperature,
    pressure,
    threshold=0.1,
    ylim=None,
    mean_percentages=None,
    reduce_compounds=True,
    init_concs=None,
    legend=True,
):
    ax = axis
    T = temperature
    P = pressure

    df = pd.DataFrame(mean_data[T][P])

    if reduce_compounds is True:
        md = df.T
        md[md["value"] != 0]
        df_p = md[md["value"] >= threshold]
        df_n = md[md["value"] <= -threshold]
        df_t = pd.concat([df_n, df_p])
        # del df_t[T][P]['CO2']
        # for x in init_concs:
        #    if not x == 'CO2':
        #        if x not in df_t[T][P] and not init_concs[x] == 0:
        #            df_t[T][P][x]['value'] = 0.0
    else:
        df_t = df

    if ylim is None:
        max_val = np.max([np.abs(df_t["value"].min()), np.abs(df_t["value"].max())])
        ylim = [-max_val - 1, max_val + 1]

    ax.set_ylim(ylim[0], ylim[1])

    norm = colors.Normalize(vmin=ylim[0], vmax=ylim[1])
    cmap = matplotlib.cm.RdBu
    colours = [cmap(norm(y)) for y in list(df_t["value"])]

    df_t.plot(kind="bar", y="value", yerr="variance", ax=ax, color=colours)

    # ax.set_yscale('symlog')
    label_names = []
    charged_species = {"CO3H": "-", "NH4": "+", "NH2CO2": "-"}
    separater = ""
    for label in list(df_t["value"].keys()):
        data = []
        for i, x in enumerate(label):
            try:
                y = int(x)
                data.append("$_{}$".format(y))
            except ValueError:
                data.append(x)
        if label in charged_species:
            data.append(str(charged_species[label]))
        label_names.append(separater.join(data))

    ax.set_ylabel("$\Delta$ Concentrations (ppm/Arb. units)")
    ax.set_xticklabels(label_names, rotation=80)

    ax.set_yticks(np.linspace(ylim[0], ylim[1], 5))
    ppm_labels = ["{:.1F}".format(x) for x in ax.get_yticks()]
    ax.set_yticklabels(ppm_labels)
    ax.axhline(0, color="k")
    rects = ax.patches
    percentages = []
    if isinstance(mean_percentages, pd.DataFrame):
        for comp, perc in mean_percentages[T][P].items():
            if not perc == np.inf:
                if math.isnan(perc) is not True:
                    percentages.append("{:.0F}".format(perc))
                else:
                    percentages.append(None)
            else:
                percentages.append(None)
    # print(percentages)
    new_label_names = []
    if not len(percentages) == 0:
        for i in range(len(label_names)):
            if percentages[i] is not None:
                new_label_names.append("{}\n{}%".format(label_names[i], percentages[i]))
            else:
                new_label_names.append(label_names[i])
    label_names = new_label_names

    for rect, label in zip(rects, label_names):
        y = rect.get_height()
        x = rect.get_x() + rect.get_width() / 2

        space = 5
        va = "bottom"
        if y < 0:
            space *= -1
            va = "top"
        if y > threshold or y < -threshold:
            ax.annotate(
                label,
                (x, y),
                xytext=(0, space),
                textcoords="offset points",
                ha="center",
                va=va,
                rotation=60,
            )
    ax.grid(axis="y", alpha=0.5)

    threshold_rectangle = patches.Rectangle(
        (0 - rects[0].get_width(), -threshold),
        len(label_names),
        threshold * 2,
        alpha=0.5,
        color="grey",
        edgecolor=None,
    )
    ax.add_patch(threshold_rectangle)
    legend_patch = patches.Patch(
        color="grey",
        alpha=0.5,
        label="detection threshold ({:.1F} ppm)".format(threshold / 1e-6),
    )
    if legend is True:
        plt.legend(handles=[legend_patch], frameon=True, edgecolor="black")
    else:
        plt.legend([], frameon=False)


def overlay_experiment(ax, exp_dict, init_concs, mean_data, T, P, threshold, **kw):
    # todo: currently doesn't plot it if included in exp_dict, but not in mean_data#
    mean_data[T][P][mean_data[T][P] != 0]
    dfp = mean_data[T][P][mean_data[T][P] >= threshold]
    dfn = mean_data[T][P][mean_data[T][P] <= -threshold]
    df_t = {T: {P: pd.concat([dfn, dfp])}}
    for x in init_concs:
        if not x == "CO2":
            if x not in df_t[T][P] and not init_concs[x] == 0:
                df_t[T][P][x] = 0.0

    available_compounds = list(df_t[T][P].keys())

    experimental_compounds = {}
    for x in available_compounds:
        if x in exp_dict:
            experimental_compounds[x] = exp_dict[x]
        else:
            experimental_compounds[x] = 0

    ax.bar(
        list(experimental_compounds), experimental_compounds.values(), fill=False, **kw
    )
    return experimental_compounds


class PrettyPlotGraph:
    def __init__(self, graph, concs, path, index, directory):
        self.graph = graph
        self.concs = concs
        self.path = path
        self.index = index
        self.directory = directory

    def make_options(self, p):
        node_sizes = [1 if isinstance(n, str) else 0 for n in self.graph.nodes()]
        node_colours = []
        alphas = []
        for n in self.graph.nodes:
            if isinstance(n, str):
                if n in list(self.concs.keys()):
                    node_colours.append((0.8, 0.0, 0.0))
                    alphas.append(1.0)
                else:
                    node_colours.append((0.6, 0.4, 0.9))
                    alphas.append(0.5)
            else:
                node_colours.append((0.0, 0.4, 0.8))
                alphas.append(0.2)

        node_options = {
            "node_color": node_colours,
            "alpha": alphas,
            "node_size": node_sizes,
        }

        edge_colours = []
        for e in self.graph.edges:
            if p[0] == e[0] and p[1] == e[1]:
                edge_colours.append((0.9, 0.0, 0.0, 1.0))
            else:
                edge_colours.append((0.2, np.random.random(), np.random.random(), 0.1))

        edge_options = {
            "connectionstyle": "arc3,rad=0.9",
            "width": 1,
            "edge_color": edge_colours,
        }

        return [node_options, edge_options]

    def make_directory(self):
        path = os.path.join(self.directory, str(self.index))
        with contextlib.suppress(FileExistsError):
            os.mkdir(path)
        return path
