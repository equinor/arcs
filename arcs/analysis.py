from monty.serialization import loadfn,dumpfn
import pandas as pd 
from collections import defaultdict
import numpy as np 
from collections import Counter
from dash import html
import warnings 
class AnalyseSampling:
    
    def __init__(
            self,
            use_markdown:bool = False,
            use_latex:bool = False,
            ):
        
        if use_markdown and use_latex:
            warnings.warn('both use_markdown and use_latex are set to True, proceeding with use_markdown=True')
        self.use_markdown = use_markdown
        self.use_latex = use_latex
            
    @staticmethod
    def latex_equation(
        equation: str,
        use_markdown: bool = True,
    ) -> str:
        r,p = equation.split('=')
        reacs = r.split(' ')
        prods = p.split(' ')
        def _latex_format(reaction_elements):
            reacs_adjusted = []
            for i in reaction_elements:
                try:
                    int(i)
                    reacs_adjusted.append(i)
                except Exception:
                    if i == '+':
                        reacs_adjusted.append(' + ')
                    else:
                        new_i = []
                        for x in i:
                            try:
                                x = int(x)
                                if use_markdown:
                                    new_i.append('<sub>{}</sub>'.format(x))
                                else:
                                    new_i.append('$_{}$'.format(x))
                            except Exception:
                                new_i.append(x)
                        reacs_adjusted.append(''.join(new_i))
            return(''.join(reacs_adjusted)) 
    
        rs = _latex_format(reacs)
        ps = _latex_format(prods)
        
        return(''.join([rs,' = ',ps]))
    
    def sci_notation(self, number, sig_fig=2):
        ret_string = "{0:.{1:d}e}".format(number, sig_fig)
        a, b = ret_string.split("e")
        # remove leading "+" and strip leading zeros
        b = int(b)
        return (a + " * 10^" + str(b))

    @staticmethod
    def get_stats(equations):

        appearances = defaultdict(int)
        for sample in equations:
            for i in sample:
                appearances[i] += 1

        equation_statistics = {}
        for equation,frequency in appearances.items():
            eq,k = equation.split(';')
            if self.cancel_markdown == True:
                equation_statistics[eq] = {'k':k.split('\n')[0],'frequency':frequency}
            else: 
                equation_statistics[self._latex_equation(eq)] = {'k':k.split('\n')[0],'frequency':frequency}
        try:
            d = pd.DataFrame(equation_statistics).T.sort_values(by='frequency',ascending=False)
            d = d.reset_index()
            d.T['index'] = 'reaction'
            d = d.to_dict()
        except Exception:
            d = {}
        return(d)
            
    def reaction_statistics(
            self,
            data:dict
            )->dict:
        """
        function that takes the reaction statistics of each sample and combines them
        returns a dict
        """
        equations = []
        for sample in data:
            for equation in sample['equation_statistics'].values():
                if equation:
                    equations.append(equation)

        statistics = Counter(equations)

        if self.use_markdown:
            statistics = {self.latex_equation(k,use_markdown=True):v for k,v in statistics.items()}
        elif self.use_latex:
            statistics = {self.latex_equation(k,use_markdown=False):v for k,v in statistics.items()}


        return(statistics)
    
    @staticmethod
    def average_sampling(
            data:dict,
    ) -> dict:
        """given a set of data, returns the initial,mean,diff and std"""
        
        average_data = {0:None,1:[]}
        for sample in data:
            for i,concentrations in sample['concentrations'].items():
                if i == 0:
                    average_data[0] = concentrations
                else:
                    average_data[1].append(concentrations)        

        import pandas as pd 
        df = pd.DataFrame(average_data[1])
        df1 = pd.DataFrame(average_data[1]).mean(axis=0)
        df0 = pd.Series(average_data[0])
        df2 = df1 - df0 
        return(pd.DataFrame([df0,df1,df2,df.std(axis=0)],index=['initial','mean','diff','std']).T)
    
    @staticmethod
    def count_path_length(
        data:dict,
    )->dict:
        return(dict(Counter([len(x['concentrations']) for x in data[1:]])))
    

    @staticmethod
    def reduce_data_by_minimum_path_length(
        data:dict,
        minimum_path_length:int = 1
    )->list:
        new_data = [data[0]]
        for _data in data[1:]:
            if len(_data['concentrations']) >=minimum_path_length:
                new_data.append(_data) 
        return(new_data)
        
        
    def reaction_paths(self,index_override=None):
        '''currently chooses the top reaction, and only does what comes after'''
        def _eqpath(pathsstats):
            _dict  = {}
            _dict['paths'] = {}
            _dict['k'] = {}
            _dict['frequency'] = pathsstats['frequency']
            for i in pathsstats['frequency']:
                r_1,k_1 = pathsstats['reaction 1'][i].split(';')
                k_1 = float(k_1.split('k=')[1])
                r_2,k_2 = pathsstats['reaction 2'][i].split(';')
                k_2 = float(k_2.split('k=')[1])
    
                str1 = r_1 + ' \n ' + r_2 
                str2 = str(k_1) + ' \n ' + str(k_2)
                _dict['paths'][i] = str1
                _dict['k'][i] = str2
            return(_dict)

        df1 = {}
        for T in self.data:
            df2 = {}
            for P in self.data[T]:
                stats = {int(x):{y:{'reaction':d.split(';')[0],'k':d.split(';')[1]} 
                                 for y,d in enumerate(self.data[T][P][x]['equation_statistics']) if d} for x in self.data[T][P]}
                                
                try:
                    if index_override == None: # should allow for clickable paths, ideally this should go through all paths
                        index = 0
                    else:
                        index = index_override
                    self.cancel_markdown=True
                    self.reaction_statistics()
                    tr = str(self.stats[float(T)][float(P)]['index'][index])
                except:
                    tr = None
    
                vs = []
                #print(stats)
                for x in stats:
                    if stats[x] and not x == 0:
                        for y in stats[x]:
                            if tr in stats[x][y]['reaction']:
                                 vs.append(x)
                                    
                self.cancel_markdown=False

                p2l = []
                for x in vs:
                    if len(stats[x]) > 1:
                        for y in stats[x]:
                            if stats[x][y]['reaction'] == tr:
                                try:
                                    r1 = self._latex_equation(stats[x][y]['reaction'])
                                    r2 = self._latex_equation(stats[x][y+1]['reaction'])
                                    #p2l.append(stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0]+':'+stats[x][y+1]['reaction']+' ; k='+stats[x][y+1]['k'].split('\n')[0])
                                    p2l.append(r1+' ; k='+stats[x][y]['k'].split('\n')[0]+':'+r2+' ; k='+stats[x][y+1]['k'].split('\n')[0])
                                except:
                                    r1 = self._latex_equation(stats[x][y-1]['reaction'])
                                    r2 = self._latex_equation(stats[x][y]['reaction'])
                                    #p2l.append(stats[x][y-1]['reaction']+' ; k='+stats[x][y-1]['k'].split('\n')[0]+':'+stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0])
                                    p2l.append(r1+' ; k='+stats[x][y-1]['k'].split('\n')[0]+':'+r2+' ; k='+stats[x][y]['k'].split('\n')[0])

                try:
                    frequencies = Counter(p2l)
                    fs = {frequencies[f]:{x:d for x,d in enumerate(f.split(':'))} for x,f in enumerate(frequencies)}
                    df = pd.DataFrame(dict(reversed(sorted(fs.items())))).T.reset_index()
    
                    df.columns = 'frequency','reaction 1','reaction 2'
                    df.set_index('frequency')
                    dict_ = df.to_dict()
                    df2[float(P)] = _eqpath(dict_)
                except:
                    df2[float(P)] = {'frequency':[None],'paths':[None],'k':[None]}

            df1[float(T)] = df2
        self.common_paths = df1
        self.cancel_markdown = False
        self.reaction_statistics()
    
    