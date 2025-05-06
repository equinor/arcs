##############old setup_functions.py
class ApplyDataToReaction:
    ''' this class applies the Gibbs data to a specific reaction'''
    
    def __init__(
            self,
            trange,
            prange,
            data,
            nprocs
            ):
        self.trange = trange
        self.prange = prange
        reactions = data['reactions']
        try:
            self.reactions = {i:Equilibrium.from_string(r) for i,r in enumerate(reactions)}
        except Exception:
            self.reactions = {i:r for i,r in enumerate(reactions)}
        self.compound_data = {k:data[k] for k in data.keys() if not k == 'reactions'}
        self.nprocs = nprocs
        self.barformat = '{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}'
        
#    def _generate_data_serial(self,t,p): #serial
#        reactions = {i:{'e':r,
#            'k':ReactionGibbsandEquilibrium(t,p,self.compound_data).#equilibrium_constant(r),
#            'g':ReactionGibbsandEquilibrium(t,p,self.compound_data).#reaction_energy(r)} 
#                     for i,r in self.reactions.items()}
#        return(reactions)
    
    def _generate_data_serial(self,t,p):
        rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
        reactions = {i:rge.as_dict(r) for i,r in self.reactions.items()}
        return(reactions)#
    def generate_data(self,t,p): #multiprocessed#
        manager = pmp.Manager()
        queue = manager.Queue()
        
        def mp_function(reaction_keys,out_q):#
            data = {}
            for r in reaction_keys:
                rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
                data[r] = rge.as_dict(self.reactions[r])
            out_q.put(data)#
        resultdict = {}
        r_keys = list(self.reactions.keys())
        chunksize = int(math.ceil(len(self.reactions)/float(self.nprocs)))
        processes = []#
        for i in range(self.nprocs):
            pr = pmp.Process(target=mp_function,
                            args=(r_keys[chunksize*i:chunksize*(i+1)],queue))
            processes.append(pr)
            pr.start()#
        for i in range(self.nprocs):
            resultdict.update(queue.get(timeout=1800))#
        for pr in processes:
            pr.join()#
        return(resultdict)#

    def apply(self,serial=False):
        data = {}
        for t in self.trange:
            pdat = {}
            for p in self.prange:
                if serial:
                    pdat[p] = self._generate_data_serial(t,p)
                else:
                    pdat[p] = self.generate_data(t,p)
            data[t] = pdat
        self.data = data
        return(self.data) 
    
    def save(self,filename='applied_reactions.p'):
        pickle.dump(self.data,open(filename,'wb'))
        print('data saved to: {}'.format(filename))




########old reactions 
import numpy as np
import itertools as it
import tqdm
from datetime import datetime
from monty.serialization import loadfn,dumpfn
import gzip,os,math
from copy import deepcopy
from pathos.helpers import mp as pmp 
from chempy import balance_stoichiometry,Reaction
from chempy.equilibria import Equilibrium
from chempy.reactionsystem import Substance
import warnings

'''
1. Starts with a ReactionDictionaryGenerator - which generates a combination of numbers
2. Filters based on compounds present
3. Filters based on balanced reactions
'''

class ReactionsDictionaryGenerator:
    ''' a class that creates the initial reference dictionary for solving all permutations of reactions between N compounds ( in this case defaults to 30 compounds) 
    
    This should initially be run at the start as it takes a long time but acts as a reference dictionary for all compounds N length or lower '''
    
    def __init__(self,no_compounds=30,path='.'):
        self.nc = no_compounds
        self.path= os.path.abspath(path)        
        
    def reaction_filter_serial(self,tl): # for serial applications
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si,sj)) or not tuple((sj,si)) in fl:
                            fl.append(tuple((si,sj)))
                    else:
                        fl.append(tuple((si,sj)))
        return(fl)
    
    def reaction_filter_mp(self,tl,out_q): # for multiprocessing 
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si,sj)) or not tuple((sj,si)) in fl:
                            fl.append(tuple((si,sj)))
                    else:
                        fl.append(tuple((si,sj)))
        out_q.put(fl)
    

    def datawriter(self,data,name):
        with open(name,'w') as f:
            [f.write(str(i)+'\n') for i in data]
            
            
    def runner_serial(self,reaction_length=4,split_files=False):
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l = self.reaction_filter_serial(tl)
                    print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l2  = self.reaction_filter_serial(tl)
                    l.extend(l2)
                    print(len(l2))
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l = self.reaction_filter_serial(tl)
                print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
    
    def _mp_run(self,tl,nprocs):
        queue = pmp.Queue()
        chunksize = int(math.ceil(len(tl)/float(nprocs)))
        
        processes = []
        for i in range(nprocs):
            pr = pmp.Process(target=self.reaction_filter_mp,
                               args=(tl[chunksize*i:chunksize*(i+1)],queue)
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
        return(data_chain)
    
    
    def _while_procs(self,tl,nprocs):
        dat = []
        while nprocs>1:
            print(nprocs,end=':')
            nprocs = int(nprocs/2)
            tl = self._mp_run(tl,nprocs)
            if len(dat) == 3:
                if dat[-2] == dat[-1] == len(tl):
                    nprocs=1
            dat.append(len(tl))
            print(len(tl),end='->')
                
        return(tl)
                
    def runner_mp(self,reaction_length=4,split_files=False,nprocs=4):
        import math
        '''protect behind if __name__ == '__main__':'''
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l = self._mp_run(tl,nprocs)
                    print(len(l))
                    #l_pre = self._mp_run(tl,nprocs)
                    #print(len(l_pre),end='->')
                    #l = self.reaction_filter_serial(l_pre)
                    #print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l2 = self._mp_run(tl,nprocs)
                    l.extend(l2)
                    print(len(l2))
                    #l_pre = self._mp_run(tl,nprocs)
                    #print(len(l_pre),end='->')
                    #l2 = self.reaction_filter_serial(l_pre)
                    #print(len(l2))
                    #l.extend(l2)
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l = self._mp_run(tl,nprocs)
                print(len(l))
               # l_pre = self._mp_run(tl,nprocs)
               #print(len(l_pre),end='->')
               # l = self.reaction_filter_serial(l_pre)
               # print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
                
    def runner_mp_while(self,reaction_length=4,split_files=False,nprocs=4):
        import math
        '''protect behind if __name__ == '__main__':'''
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')                    
                    l_pre = self._while_procs(tl,nprocs)
                    l = self.reaction_filter_serial(l_pre)
                    print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')                    
                    l_pre = self._while_procs(tl,nprocs)
                    l2 = self.reaction_filter_serial(l_pre)
                    print(len(l2))
                    l.extend(l2)
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l_pre = self._while_procs(tl,nprocs)
                l = self.reaction_filter_serial(l_pre)
                print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
                    

    def clean_from_file(self,filename,nprocs=4):
        import gzip
        '''file is a gzip'''
        def _convert_line(file):
            for line in file:
                r,p = line.strip().split('],')
                r = tuple(int(x) for x in r.split('([')[1].split(',') if x)
                p = tuple(int(x) for x in p.split('[')[1].split('])')[0].split(',') if x)
                yield ((r,p))
        print('reading from file...',end=' ')
        f = tuple(_convert_line(gzip.open(filename,'rt')))
        print(len(f),end='->')
        l = self._mp_run(f,nprocs)
        print(len(l))
        filename = filename.replace('.dat.gz','_reloaded.dat')
        self.datawriter(data=l,name=filename)




class MappingtoReaction:
    ''' a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered'''
    
    def __init__(self,filename,compounds):
        self.filename = filename
        self.compounds = {i:c for i,c in enumerate(compounds)}
        
    def convert_file(self,file): # iterator that reads in the file (file is a gzip)
        for line in file:
            r,p = line.strip().split('],')
            r = tuple(int(x) for x in r.split('([')[1].split(',') if x)
            p = tuple(int(x) for x in p.split('[')[1].split('])')[0].split(',') if x)    
            yield ((r,p))
            
    def remove_indices(self,reaction_indexes):
        ''''''
        indexes =  tuple(self.compounds)
        approved = []
        for r in reaction_indexes:
            if not [x for x in it.chain(*r) if not x in indexes]:
                approved.append(r)
        return(approved)
    
    def convert_to_string(self,approved_list):
    
        def _screen_string(r):
            re,pr = r
            c = []
            for i in re:
                rs = [x for x in i]
                for j in rs:
                    try:
                        int(j)
                    except Exception:
                        c.append(j)
        
            dr = set(dict.fromkeys(c))
        
            c = []
            for i in pr:
                ps = [x for x in i]
                for j in ps:
                    try:
                        int(j)
                    except Exception:
                        c.append(j)
        
            dp = set(dict.fromkeys(c))
        
            if dr == dp:
                return(r)
    
        converted = []
        for r in approved_list:
            d = [[self.compounds[i] for i in r[0]],[self.compounds[i] for i in r[1]]]
            ds = _screen_string(d)
            if ds:
                converted.append(ds)
        return(converted)
    
    def convert_to_equation(self,converted_strings):
        converted = []
        for r in tqdm.tqdm(converted_strings):
            re,pr = r
            try:
                converted.append(balance_stoichiometry(list(re),list(pr),underdetermined=None))
            except Exception:
                try:
                    converted.append(balance_stoichiometry(list(re),list(pr)))
                except Exception:
                    pass                            
        return(converted) 
    
    def screen_converted(self,converted_reactions): # this should be multiprocessed - > perhaps a mp.Pool? as it only needs the big list
        def _convert_ord_to_dict(r):
            re,pr =  r
            try:
                reacs = {k:int(re[k]) for k in re}
                prods = {k:int(pr[k]) for k in pr}
            except Exception:
                warnings.warn('\n error with {}'.format(r))
            return(reacs,prods)
    
        screened = []
        for r in converted_reactions:
            try:
                re,pr = _convert_ord_to_dict(r)
                try:
                    screened.append(Equilibrium(re,pr))
                except Exception:
                    pass
            except Exception :
                pass
        return(screened)    
    
    def run_all(self):
        s = datetime.now()
        loi = tuple(self.convert_file(gzip.open(self.filename,'rt')))
        print('orig =',len(loi),end='...')
        approved = self.remove_indices(loi)
        print(' approved = ',len(approved),end='...')
        strings = self.convert_to_string(approved)
        print(' prescreening = ',len(strings))
        equations = self.convert_to_equation(strings)
        screened = self.screen_converted(equations)
        print(' final = ',len(screened),end='...')
        f = datetime.now() - s 
        print(' time = ',f)
        return(screened)
    


################################traversal


    def _random_choice_unconnected(self,T,P,force_direct=False,co2=False): # currently randomly disjointed reactions that are weighted
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        if force_direct:
            pstring = [0,1,2]
            while len(pstring) > 2:
                source = self._get_weighted_random_compound(T,P,co2=co2,force_selection=None) 
                target = np.random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target,weight='weight')
                pstring = [n for n in p if isinstance(p,str)]
        else:
                source = self._get_weighted_random_compound(T,P)
                target = np.random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
        return(p)

    def _random_choice_connected(self,T,P,force_direct=False,previous_index=None,co2=False): # this will be joined - I think we can also make a ranking of potential reactions based upon components in the stream as well 
        if previous_index == None:
            raise ValueError('no previous compound selected')
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        if force_direct:
            pstring = [0,1,2]
            while len(pstring) > 2:
                present = [c for c in list(self.reactions[T][P][previous_index]['e'].reac) + list(self.reactions[T][P][previous_index]['e'].prod) ] # this should probably be weighted according to stoichiometry i.e. 2CO2 + H2O = [CO2, CO2, H2O]
                source = self._get_weighted_random_compound(T,P,co2=co2,force_selection=present)
                target = np.random.choice(nodes) # the next path will be random 
                p = nx.shortest_path(self.graph[T][P],source,target,weight='weight')
                pstring = [n for n in p if isinstance(p,str)]
        else:
                source = self._get_weighted_random_compound(T,P)
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
        return(p)
