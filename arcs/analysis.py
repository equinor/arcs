from monty.serialization import loadfn
import pandas as pd
from collections import defaultdict
import numpy as np
from collections import Counter


class AnalyseSampling:
    def __init__(self, data, markdown=False):
        if isinstance(data, str):
            self.data = loadfn(data)
        else:
            self.data = data
        self.markdown = markdown

    def _latex_equation(self, equation):
        r, p = equation.split("=")
        reacs = r.split(" ")
        prods = p.split(" ")

        def _latex_format(reaction_elements):
            reacs_adjusted = []
            for i in reaction_elements:
                try:
                    int(i)
                    reacs_adjusted.append(i)
                except ValueError:
                    if i == "+":
                        reacs_adjusted.append(" + ")
                    else:
                        new_i = []
                        for x in i:
                            try:
                                x = int(x)
                                if self.markdown is True:
                                    new_i.append("<sub>{}</sub>".format(x))
                                else:
                                    new_i.append("$_{}$".format(x))
                            except ValueError:
                                new_i.append(x)
                        reacs_adjusted.append("".join(new_i))
            return "".join(reacs_adjusted)

        rs = _latex_format(reacs)
        ps = _latex_format(prods)
        return "".join([rs, " = ", ps])

    def _sci_notation(self, number, sig_fig=2):
        ret_string = "{0:.{1:d}e}".format(number, sig_fig)
        a, b = ret_string.split("e")
        # remove leading "+" and strip leading zeros
        b = int(b)
        return a + " * 10^" + str(b)

    def _get_stats(self, equations):
        appearances = defaultdict(int)
        for sample in equations:
            for i in sample:
                appearances[i] += 1

        equation_statistics = {}
        for equation, frequency in appearances.items():
            eq, k = equation.split(";")
            if self.cancel_markdown is True:
                equation_statistics[eq] = {
                    "k": k.split("\n")[0],
                    "frequency": frequency,
                }
            else:
                equation_statistics[self._latex_equation(eq)] = {
                    "k": k.split("\n")[0],
                    "frequency": frequency,
                }
        try:
            d = pd.DataFrame(equation_statistics).T.sort_values(
                by="frequency", ascending=False
            )
            d = d.reset_index()
            d.T["index"] = "reaction"
            d = d.to_dict()
        except Exception:
            d = {}
        return d

    def reaction_statistics(self):
        eqt = {}
        for T in self.data:
            eqp = {}
            for P in self.data[T]:
                equations = []
                for x in self.data[T][P]:
                    eqs = self.data[T][P][x]["equation_statistics"]
                    if eqs:
                        equations.append(eqs)
                eqp[P] = self._get_stats(equations)
            eqt[T] = eqp

        self.stats = eqt

    def mean_sampling(self):
        t = list(self.data)[0]
        p = list(self.data[t])[0]
        zeroth = list(self.data[t][p])[0]
        # if isinstance(xr[0],str):
        #    #data = str_to_int_dict(data)

        md = {}
        f_1 = {}

        for T in self.data:
            f_2 = {}
            m_2 = {}
            for P in self.data[T]:
                f_3 = {}
                m_3 = {}
                for c in self.data[T][P][zeroth]["data"]:
                    f_4 = []
                    for x in self.data[T][P]:
                        if not x == zeroth:
                            f_4.append(self.data[T][P][x]["data"][c])
                    diff = [i - self.data[T][P][zeroth]["data"][c] for i in f_4]
                    f_3[c] = np.mean(f_4) / 1e-6
                    m_3[c] = {
                        "value": np.mean(diff) / 1e-6,
                        "variance": np.var(diff) / 1e-6,
                    }  # 2nd value is the variance and not the std deviation
                m_2[P] = m_3
                f_2[P] = f_3
            md[T] = m_2
            f_1[T] = f_2

        self.mean_data = md
        self.final_concs = f_1

    def reaction_paths(self, index_override=None):
        """currently chooses the top reaction, and only does what comes after"""

        def _eqpath(pathsstats):
            _dict = {}
            _dict["paths"] = {}
            _dict["k"] = {}
            _dict["frequency"] = pathsstats["frequency"]
            for i in pathsstats["frequency"]:
                r_1, k_1 = pathsstats["reaction 1"][i].split(";")
                k_1 = float(k_1.split("k=")[1])
                r_2, k_2 = pathsstats["reaction 2"][i].split(";")
                k_2 = float(k_2.split("k=")[1])

                str1 = r_1 + " \n " + r_2
                str2 = str(k_1) + " \n " + str(k_2)
                _dict["paths"][i] = str1
                _dict["k"][i] = str2
            return _dict

        df1 = {}
        for T in self.data:
            df2 = {}
            for P in self.data[T]:
                stats = {
                    int(x): {
                        y: {"reaction": d.split(";")[0], "k": d.split(";")[1]}
                        for y, d in enumerate(self.data[T][P][x]["equation_statistics"])
                        if d
                    }
                    for x in self.data[T][P]
                }

                try:
                    if (
                        index_override is None
                    ):  # should allow for clickable paths, ideally this should go through all paths
                        index = 0
                    else:
                        index = index_override
                    self.cancel_markdown = True
                    self.reaction_statistics()
                    tr = str(self.stats[T][P]["index"][index])
                except Exception:
                    tr = None

                vs = []
                # print(stats)
                for x in stats:
                    if stats[x] and not x == 0:
                        for y in stats[x]:
                            if tr in stats[x][y]["reaction"]:
                                vs.append(x)

                self.cancel_markdown = False

                p2l = []
                for x in vs:
                    if len(stats[x]) > 1:
                        for y in stats[x]:
                            if stats[x][y]["reaction"] == tr:
                                try:
                                    r1 = self._latex_equation(stats[x][y]["reaction"])
                                    r2 = self._latex_equation(
                                        stats[x][y + 1]["reaction"]
                                    )
                                    # p2l.append(stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0]+':'+stats[x][y+1]['reaction']+' ; k='+stats[x][y+1]['k'].split('\n')[0])
                                    p2l.append(
                                        r1
                                        + " ; k="
                                        + stats[x][y]["k"].split("\n")[0]
                                        + ":"
                                        + r2
                                        + " ; k="
                                        + stats[x][y + 1]["k"].split("\n")[0]
                                    )
                                except Exception:
                                    r1 = self._latex_equation(
                                        stats[x][y - 1]["reaction"]
                                    )
                                    r2 = self._latex_equation(stats[x][y]["reaction"])
                                    # p2l.append(stats[x][y-1]['reaction']+' ; k='+stats[x][y-1]['k'].split('\n')[0]+':'+stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0])
                                    p2l.append(
                                        r1
                                        + " ; k="
                                        + stats[x][y - 1]["k"].split("\n")[0]
                                        + ":"
                                        + r2
                                        + " ; k="
                                        + stats[x][y]["k"].split("\n")[0]
                                    )

                try:
                    frequencies = Counter(p2l)
                    fs = {
                        frequencies[f]: {x: d for x, d in enumerate(f.split(":"))}
                        for x, f in enumerate(frequencies)
                    }
                    df = pd.DataFrame(
                        dict(reversed(sorted(fs.items())))
                    ).T.reset_index()

                    df.columns = "frequency", "reaction 1", "reaction 2"
                    df.set_index("frequency")
                    dict_ = df.to_dict()
                    df2[P] = _eqpath(dict_)
                except Exception:
                    df2[P] = {"frequency": [None], "paths": [None], "k": [None]}

            df1[T] = df2
        self.common_paths = df1
        self.cancel_markdown = False
        self.reaction_statistics()
