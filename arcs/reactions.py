import numpy as np
import itertools as it
import tqdm
from datetime import datetime
import gzip
import os
import math
from pathos.helpers import mp as pmp
from chempy import balance_stoichiometry
from chempy.equilibria import Equilibrium
import warnings

"""
1. Starts with a ReactionDictionaryGenerator - which generates a combination of numbers
2. Filters based on compounds present
3. Filters based on balanced reactions
"""


class ReactionsDictionaryGenerator:
    """a class that creates the initial reference dictionary for solving all permutations of reactions between N compounds ( in this case defaults to 30 compounds)

    This should initially be run at the start as it takes a long time but acts as a reference dictionary for all compounds N length or lower"""

    def __init__(self, no_compounds=30, path="."):
        self.nc = no_compounds
        self.path = os.path.abspath(path)

    def reaction_filter_serial(self, tl):  # for serial applications
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si, sj)) or tuple((sj, si)) not in fl:
                            fl.append(tuple((si, sj)))
                    else:
                        fl.append(tuple((si, sj)))
        return fl

    def reaction_filter_mp(self, tl, out_q):  # for multiprocessing
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si, sj)) or tuple((sj, si)) not in fl:
                            fl.append(tuple((si, sj)))
                    else:
                        fl.append(tuple((si, sj)))
        out_q.put(fl)

    def datawriter(self, data, name):
        with open(name, "w") as f:
            [f.write(str(i) + "\n") for i in data]

    def runner_serial(self, reaction_length=4, split_files=False):
        sizing = [
            x
            for x in it.combinations_with_replacement(np.arange(1, reaction_length), 2)
            if np.sum(x) == reaction_length
        ]

        print(
            "running_runner for {} component reactions...".format(reaction_length),
            end=" ",
        )
        print("{} possibilities : {}".format(len(sizing), sizing))

        if len(sizing) == 1 or split_files is False:
            for i, size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    data = self.reaction_filter_serial(tl)
                    print(len(data))
                else:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    l2 = self.reaction_filter_serial(tl)
                    data.extend(l2)
                    print(len(l2))

                filename = os.path.join(
                    self.path, "data_{}.dat".format(reaction_length)
                )
                self.datawriter(data=data, name=filename)
        else:
            for i, size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                tl = tuple(it.product(rs, ps))
                print(size, ":", len(tl), end="->")
                data = self.reaction_filter_serial(tl)
                print(len(data))
                filename = os.path.join(
                    self.path, "data_{}-{}.dat".format(reaction_length, i)
                )
                self.datawriter(data=data, name=filename)

    def _mp_run(self, tl, nprocs):
        queue = pmp.Queue()
        chunksize = int(math.ceil(len(tl) / float(nprocs)))

        processes = []
        for i in range(nprocs):
            pr = pmp.Process(
                target=self.reaction_filter_mp,
                args=(tl[chunksize * i : chunksize * (i + 1)], queue),
            )
            processes.append(pr)
            pr.start()

        data = []
        for i in range(nprocs):
            data.append(queue.get())

        for pr in processes:
            pr.join()

        data_chain = tuple(it.chain(*data))
        # should the serial part come here?
        return data_chain

    def _while_procs(self, tl, nprocs):
        dat = []
        while nprocs > 1:
            print(nprocs, end=":")
            nprocs = int(nprocs / 2)
            tl = self._mp_run(tl, nprocs)
            if len(dat) == 3:
                if dat[-2] == dat[-1] == len(tl):
                    nprocs = 1
            dat.append(len(tl))
            print(len(tl), end="->")

        return tl

    def runner_mp(self, reaction_length=4, split_files=False, nprocs=4):
        """protect behind if __name__ == '__main__':"""
        sizing = [
            x
            for x in it.combinations_with_replacement(np.arange(1, reaction_length), 2)
            if np.sum(x) == reaction_length
        ]

        print(
            "running_runner for {} component reactions...".format(reaction_length),
            end=" ",
        )
        print("{} possibilities : {}".format(len(sizing), sizing))

        if len(sizing) == 1 or split_files is False:
            for i, size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    data = self._mp_run(tl, nprocs)
                    print(len(data))
                    # l_pre = self._mp_run(tl,nprocs)
                    # print(len(l_pre),end='->')
                    # l = self.reaction_filter_serial(l_pre)
                    # print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    l2 = self._mp_run(tl, nprocs)
                    data.extend(l2)
                    print(len(l2))
                    # l_pre = self._mp_run(tl,nprocs)
                    # print(len(l_pre),end='->')
                    # l2 = self.reaction_filter_serial(l_pre)
                    # print(len(l2))
                    # l.extend(l2)

                filename = os.path.join(
                    self.path, "data_{}.dat".format(reaction_length)
                )
                self.datawriter(data=data, name=filename)
        else:
            for i, size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                tl = tuple(it.product(rs, ps))
                print(size, ":", len(tl), end="->")
                data = self._mp_run(tl, nprocs)
                print(len(data))
                # l_pre = self._mp_run(tl,nprocs)
                # print(len(l_pre),end='->')
                # l = self.reaction_filter_serial(l_pre)
                # print(len(l))
                filename = os.path.join(
                    self.path, "data_{}-{}.dat".format(reaction_length, i)
                )
                self.datawriter(data=data, name=filename)

    def runner_mp_while(self, reaction_length=4, split_files=False, nprocs=4):
        """protect behind if __name__ == '__main__':"""
        sizing = [
            x
            for x in it.combinations_with_replacement(np.arange(1, reaction_length), 2)
            if np.sum(x) == reaction_length
        ]

        print(
            "running_runner for {} component reactions...".format(reaction_length),
            end=" ",
        )
        print("{} possibilities : {}".format(len(sizing), sizing))

        if len(sizing) == 1 or split_files is False:
            for i, size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    l_pre = self._while_procs(tl, nprocs)
                    data = self.reaction_filter_serial(l_pre)
                    print(len(data))
                else:
                    rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                    ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                    tl = tuple(it.product(rs, ps))
                    print(size, ":", len(tl), end="->")
                    l_pre = self._while_procs(tl, nprocs)
                    l2 = self.reaction_filter_serial(l_pre)
                    print(len(l2))
                    data.extend(l2)

                filename = os.path.join(
                    self.path, "data_{}.dat".format(reaction_length)
                )
                self.datawriter(data=data, name=filename)
        else:
            for i, size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc + 1)], size[0])
                ps = it.combinations([x for x in range(self.nc + 1)], size[1])
                tl = tuple(it.product(rs, ps))
                print(size, ":", len(tl), end="->")
                l_pre = self._while_procs(tl, nprocs)
                data = self.reaction_filter_serial(l_pre)
                print(len(data))
                filename = os.path.join(
                    self.path, "data_{}-{}.dat".format(reaction_length, i)
                )
                self.datawriter(data=data, name=filename)

    def clean_from_file(self, filename, nprocs=4):
        import gzip

        """file is a gzip"""

        def _convert_line(file):
            for line in file:
                r, p = line.strip().split("],")
                r = tuple(int(x) for x in r.split("([")[1].split(",") if x)
                p = tuple(
                    int(x) for x in p.split("[")[1].split("])")[0].split(",") if x
                )
                yield ((r, p))

        print("reading from file...", end=" ")
        f = tuple(_convert_line(gzip.open(filename, "rt")))
        print(len(f), end="->")
        data_chain = self._mp_run(f, nprocs)
        print(len(data_chain))
        filename = filename.replace(".dat.gz", "_reloaded.dat")
        self.datawriter(data=data_chain, name=filename)


class MappingtoReaction:
    """a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered"""

    def __init__(self, filename, compounds):
        self.filename = filename
        self.compounds = {i: c for i, c in enumerate(compounds)}

    def convert_file(self, file):  # iterator that reads in the file (file is a gzip)
        for line in file:
            r, p = line.strip().split("],")
            r = tuple(int(x) for x in r.split("([")[1].split(",") if x)
            p = tuple(int(x) for x in p.split("[")[1].split("])")[0].split(",") if x)
            yield ((r, p))

    def remove_indices(self, reaction_indexes):
        indexes = tuple(self.compounds)
        approved = []
        for r in reaction_indexes:
            if not [x for x in it.chain(*r) if x not in indexes]:
                approved.append(r)
        return approved

    def convert_to_string(self, approved_list):
        def _screen_string(r):
            re, pr = r
            c = []
            for i in re:
                rs = [x for x in i]
                for j in rs:
                    try:
                        int(j)
                    except ValueError:
                        c.append(j)

            dr = set(dict.fromkeys(c))

            c = []
            for i in pr:
                ps = [x for x in i]
                for j in ps:
                    try:
                        int(j)
                    except ValueError:
                        c.append(j)

            dp = set(dict.fromkeys(c))

            if dr == dp:
                return r

        converted = []
        for r in approved_list:
            d = [[self.compounds[i] for i in r[0]], [self.compounds[i] for i in r[1]]]
            ds = _screen_string(d)
            if ds:
                converted.append(ds)
        return converted

    def convert_to_equation(self, converted_strings):
        converted = []
        for r in tqdm.tqdm(converted_strings):
            re, pr = r
            try:
                converted.append(
                    balance_stoichiometry(list(re), list(pr), underdetermined=None)
                )
            except Exception:
                try:
                    converted.append(balance_stoichiometry(list(re), list(pr)))
                except Exception:
                    pass
        return converted

    def screen_converted(
        self, converted_reactions
    ):  # this should be multiprocessed - > perhaps a mp.Pool? as it only needs the big list
        def _convert_ord_to_dict(r):
            re, pr = r
            try:
                reacs = {k: int(re[k]) for k in re}
                prods = {k: int(pr[k]) for k in pr}
            except Exception:
                warnings.warn("\n error with {}".format(r))
            return (reacs, prods)

        screened = []
        for r in converted_reactions:
            try:
                re, pr = _convert_ord_to_dict(r)
                try:
                    screened.append(Equilibrium(re, pr))
                except Exception:
                    pass
            except Exception:
                pass
        return screened

    def run_all(self):
        s = datetime.now()
        loi = tuple(self.convert_file(gzip.open(self.filename, "rt")))
        print("orig =", len(loi), end="...")
        approved = self.remove_indices(loi)
        print(" approved = ", len(approved), end="...")
        strings = self.convert_to_string(approved)
        print(" prescreening = ", len(strings))
        equations = self.convert_to_equation(strings)
        screened = self.screen_converted(equations)
        print(" final = ", len(screened), end="...")
        f = datetime.now() - s
        print(" time = ", f)
        return screened
