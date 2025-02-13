from typing import Any, Hashable, no_type_check
import pandas as pd
from collections import defaultdict
import numpy as np
from collections import Counter

from .traversal import _RandomWalk


class AnalyseSampling:
    def __init__(self, data: dict[int, _RandomWalk]) -> None:
        self.data = data

    def _latex_equation(self, equation: str) -> str:
        r, p = equation.split("=")
        reacs = r.split(" ")
        prods = p.split(" ")

        def _latex_format(reaction_elements: list[str]) -> str:
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
                                new_i.append("<sub>{}</sub>".format(int(x)))
                            except ValueError:
                                new_i.append(x)
                        reacs_adjusted.append("".join(new_i))
            return "".join(reacs_adjusted)

        rs = _latex_format(reacs)
        ps = _latex_format(prods)
        return "".join([rs, " = ", ps])

    def _sci_notation(self, number: float, sig_fig: int = 2) -> str:
        ret_string = "{0:.{1:d}e}".format(number, sig_fig)
        a, b = ret_string.split("e")
        # remove leading "+" and strip leading zeros
        return a + " * 10^" + str(int(b))

    def _get_stats(self, equations: list[list[str]]) -> dict[Hashable, Any]:
        appearances: dict[str, int] = defaultdict(int)
        for sample in equations:
            for i in sample:
                appearances[i] += 1

        equation_statistics = {}
        for equation, frequency in appearances.items():
            eq, k = equation.split(";")
            equation_statistics[eq] = {
                "k": k.split("\n")[0],
                "frequency": frequency,
            }
        try:
            d = pd.DataFrame(equation_statistics).T.sort_values(
                by="frequency", ascending=False
            )
            d = d.reset_index()
            d.T["index"] = "reaction"
            return d.to_dict()
        except Exception:
            pass
        return {}

    def reaction_statistics(self) -> None:
        equations: list[list[str]] = []
        for x in self.data:
            eqs = self.data[x]["equation_statistics"]
            if eqs:
                equations.append(eqs)
        self.stats = self._get_stats(equations)

    def mean_sampling(self) -> None:
        final_concentrations = {}
        mean_values = {}

        compounds: set[str] = set()
        for sample in self.data.values():
            compounds.update(sample["data"].keys())

        # initial concentrations including zeros
        initial_concentrations = {c: self.data[0]["data"].get(c, 0) for c in compounds}

        for compound in compounds:
            final_all_samples = [
                self.data[x]["data"].get(compound, 0) for x in self.data
            ]
            diff = [i - initial_concentrations[compound] for i in final_all_samples]
            final_concentrations[compound] = np.mean(final_all_samples) / 1e-6
            mean_values[compound] = {
                "value": np.mean(diff) / 1e-6,
                "variance": np.var(diff) / 1e-6,
            }
            # 2nd value is the variance and not the std deviation

        self.mean_data = mean_values
        self.final_concs = final_concentrations

    def reaction_paths(self, index_override: int | None = None) -> None:
        """currently chooses the top reaction, and only does what comes after"""

        @no_type_check
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

        stats = {
            int(x): {
                y: {"reaction": d.split(";")[0], "k": d.split(";")[1]}
                for y, d in enumerate(self.data[x]["equation_statistics"])
                if d
            }
            for x in self.data
        }

        try:
            if (
                index_override is None
            ):  # should allow for clickable paths, ideally this should go through all paths
                index = 0
            else:
                index = index_override
            self.reaction_statistics()
            tr = str(self.stats["index"][index])
        except Exception:
            tr = None

        vs = []
        # print(stats)
        for x in stats:
            if stats[x] and not x == 0:
                for y in stats[x]:
                    if tr in stats[x][y]["reaction"]:
                        vs.append(x)

        p2l = []
        for x in vs:
            if len(stats[x]) > 1:
                for y in stats[x]:
                    if stats[x][y]["reaction"] == tr:
                        try:
                            r1 = self._latex_equation(stats[x][y]["reaction"])
                            r2 = self._latex_equation(stats[x][y + 1]["reaction"])
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
                            r1 = self._latex_equation(stats[x][y - 1]["reaction"])
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
            df = pd.DataFrame(dict(reversed(sorted(fs.items())))).T.reset_index()

            df.columns = pd.Index(("frequency", "reaction 1", "reaction 2"))
            df.set_index("frequency")
            dict_ = df.to_dict()
            df1 = _eqpath(dict_)
        except Exception:
            df1 = {"frequency": None, "paths": None, "k": None}

        self.common_paths = df1
        self.reaction_statistics()
