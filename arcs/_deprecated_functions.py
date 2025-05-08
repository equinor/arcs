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
    
   def _queue_function(self,
                        pbari,
                        samples,T,P,
                        probability_threshold,
                        path_depth,
                        max_compounds,
                        max_rank,
                        co2,
                        scale_highest,
                        ceiling,
                        method,
                        out_q):
        
        sample_data = {}
        with tqdm(total=len(samples),bar_format='progress: {desc:<10}|{bar:50}|',ascii=' >=',position=0,leave=False) as pbar:
            for sample in samples:
                sample_data[sample] = self.random_walk(T=T,P=P,
                                                       probability_threshold=probability_threshold,
                                                       path_depth=path_depth,
                                                       max_compounds=max_compounds,
                                                       max_rank=max_rank,
                                                       co2=co2,
                                                       scale_highest=scale_highest,
                                                       ceiling=ceiling,
                                                       method=method)
                pbar.update(1)
                    
        out_q.put(sample_data)

    
    def sampling_multiprocessing(self,T=None,P=None,**kw):
        
        init_concs = copy.deepcopy(self.concs)
        result_dict = {0:{'data':init_concs,'equation_statistics':[],'path_length':None}}
        #start the queue
        out_queue = pmp.Queue()
        samples = list(range(1,self.sample_length+1,1))
        data_chunks = [samples[chunksize*i:chunksize*(i+1)] 
                            for i in range(self.nprocs) 
                            for chunksize in [int(math.ceil(len(samples)/float(self.nprocs)))]]
        
        jobs = []
        for i,chunk in enumerate(data_chunks):
            process = pmp.Process(target=self._queue_function,
                                  args=(i,chunk,T,P,
                                        self.probability_threshold,
                                        self.path_depth,
                                        self.max_compounds,
                                        self.max_rank,
                                        self.co2,
                                        self.scale_highest,
                                        self.ceiling,
                                        self.method,
                                        out_queue))
            jobs.append(process)
            process.start()


        for proc in jobs:
            result_dict.update(out_queue.get())
            
        for proc in jobs:
            proc.terminate()

        for proc in jobs:
            proc.join()

        out_queue.close()
        


        return(result_dict) 
    
    def sampling_serial(self,T=None,P=None,**kw):
        init_concs = copy.deepcopy(self.concs)
        result_dict = {0:{'data':init_concs,'equation_statistics':[],'path_length':None}}
        with tqdm(total=self.sample_length,bar_format='progress: {desc:<10}|{bar:50}|',ascii=' >=',position=0,leave=False) as pbar:
            for sample in range(self.sample_length):
                result_dict[sample+1] = self.random_walk(T=T,P=P,
                                                       probability_threshold=self.probability_threshold,
                                                       path_depth=self.path_depth,
                                                       max_compounds=self.max_compounds,
                                                       max_rank=self.max_rank,
                                                       co2=self.co2,
                                                       scale_highest=self.scale_highest,
                                                       ceiling=self.ceiling,
                                                       method=self.method)
                pbar.update(1)
        return(result_dict)


    def run(self,trange,prange,ic=None,save=False,savename=None,ignore_warnings=True,logging=False,**kw):
        if ignore_warnings==True:
            warnings.filterwarnings("ignore")

        '''
        kwargs = sample_length,probability_threshold,max_compounds,max_rank,path_depth,nprocs,random_path_depth,co2=False
        '''
        from loguru import logger
        from io import StringIO
        #setup logger
        stream = StringIO()
        logger.remove()
        logger.add(stream,format="{message}")

        num=1
        total = len(trange) * len(prange)
        
        from datetime import datetime
        needed_args = self.__dict__
        for i in needed_args:
            if i in kw:
                self.__dict__[i] = kw[i]
            
        logger.info('''\n                                             
                                            
    // | |     //   ) )  //   ) )  //   ) ) 
   //__| |    //___/ /  //        ((        
  / ___  |   / ___ (   //           \\      
 //    | |  //   | |  //              ) )   
//     | | //    | | ((____/ / ((___ / /    
version:1.2
{}
        ->sample_length = {}
        ->probability_threshold = {}
        ->max_compounds = {}
        ->max_rank = {}
        ->path_depth = {}
        ->co2 = {}
        ->shortest path method = {}
        ->number of processes = {}
        ->concentration ceiling = {} %
        ->scale highest = {}
        ->rank smaller reactions higher = {}\n'''.format(str(datetime.now()),self.sample_length,
                                       self.probability_threshold,self.max_compounds,
                                       self.max_rank,self.path_depth,self.co2,self.method,
                                       self.nprocs,self.ceiling,self.scale_highest,self.rank_small_reactions_higher))
        
        logger.info('initial concentrations (ppm):\n')
        self.concs = ic
        concstring = pd.Series({k:v for k,v, in self.concs.items() if v > 0}) / 1e-6
        del concstring['CO2']
        logger.info(concstring.to_string()+'\n')

        if logging:
            print(stream.get_value())
        
        
        path_lengths = [] 
        total_data = {}
        for T in trange:
            data_2 = {}
            final_concs_2 = {}
            initfinaldiff = {}
            for P in prange:
                start = datetime.now()
                logger.info('\n {}/{}: temperature = {}K, pressure = {}bar '.format(num,total,T,P),end='\n')
                if self.nprocs > 1:
                    data_2[P] =  self.sampling_multiprocessing(T,P,**kw)
                else:
                    data_2[P] = self.sampling_serial(T,P,**kw)

                finish = datetime.now() - start
                logger.info('-> completed in {} seconds'.format(finish.total_seconds()),end='\n')
                reformatted = [{x:v for x,v in data_2[P][i]['data'].items()} for i in data_2[P]]
                mean = pd.Series({k:v for k,v in pd.DataFrame(reformatted).mean().items() if v > 0.5e-6}).drop('CO2')/1e-6
                #mean = pd.Series({x:v for x,v in np.mean(pd.DataFrame(data_2[P][i]['data'] for i in data_2[P]).keys()) if v > 0.5e-6}).drop('CO2')/1e-6
                logger.info('\n final concentrations (>0.5ppm):\n')
                logger.info(mean.round(1).to_string())
                final_concs_2[P] = mean.to_dict()
                diff_concs = pd.Series(mean.to_dict()) - pd.Series({k:v/1e-6 for k,v in self.concs.items()})
                ift = pd.DataFrame([{k:v/1e-6 for k,v in self.concs.items() if v > 0},mean.to_dict(),diff_concs.to_dict()],index=['initial','final','change']).T
                initfinaldiff[P] = ift.dropna(how='all').fillna(0.0).to_dict()
                avgpathlength = np.median([data_2[P][i]['path_length'] for i in data_2[P] if not data_2[P][i]['path_length'] == None])


                logger.info('\n median path length: {}'.format(avgpathlength))
                path_lengths.append(avgpathlength)
                num+=1
            total_data[T] = data_2
            self.final_concs[T] = final_concs_2
            self.initfinaldiff[T] = initfinaldiff            
                
        if save:
            from monty.serialization import dumpfn
            if not savename:
                from datetime import date
                today = str(date.today())
                savename='sampling_{}.json'.format(today)
            dumpfn(total_data,savename,indent=4)

        self.metadata = {'arcs_version':"1.4.0",
                         'avg_path_length':np.mean(path_lengths),
                         'co2':self.co2,
                         'max_compounds':self.max_compounds,
                         'probability_threshold':self.probability_threshold,
                         'shortest_path_method':self.method,
                         'max_rank':self.max_rank,
                         'sample_length':self.sample_length,
                         'path_depth':self.path_depth,
                         'random_path_depth':self.random_path_depth,
                         'nprocs':self.nprocs,
                         'ceiling':self.ceiling,
                         'scale_highest':self.scale_highest,
                         'rank_small_reactions_higher':self.rank_small_reactions_higher,
                         'platform':platform.platform(),
                         'python_version':platform.python_version(),
                         'processor':platform.processor(),
                         'available_cores':psutil.cpu_count(),
                         'available_memory':str(int(psutil.virtual_memory()[0] / 1000/1000/1000))+'Gb',
                         'date':str(datetime.now())}
        
       
        self.data = total_data  
        if logging:
            print(stream.get_value())    

    def _get_weighted_reaction_rankings(
            self,
            weighted_random_compounds:dict,
            maximum_reaction_number:int = 20,
            shortest_path_method:str = 'Djikstra',
            ):
        """
        given a dictionary of weighted random compounds from self.get_weighted_random_compounds find the shortest path between: 
        1. compound 0 and compound 1 -> list
        2. filter based on available compounds 
        3. weight based on edge weight and (optional) length multiplier for unreasonable large reactions (if option selected)
        returns a dictionary of ranked reactions and their weighting. 
        """
        
        if len(weighted_random_compounds) == 1:
            return(None)

        rankings = {}

        shortest_paths = nx.shortest_paths.all_shortest_paths(
                G=self.graph,
                source=list(weighted_random_compounds)[0],
                target=list(weighted_random_compounds)[1],
                method=shortest_path_method)  # gets a list of shortest paths without weights first
                
        rankings = {}
        for path in shortest_paths:
            source,reaction,target = path
            reaction_compounds = list(self.graph[reaction])
            #find reactions with >3rd compound
            if len(weighted_random_compounds) > 2:
                for compound in list(weighted_random_compounds)[2:]:
                    if compound in reaction_compounds:
                        rankings[reaction] = self.graph.get_edge_data(
                            u=source,
                            v=reaction
                        )[0]['weight']*self.length_multiplier(reaction) #need to play around with coefficients=True
            else:
                rankings[reaction] = self.graph.get_edge_data(
                    u=source,
                    v=reaction
                )[0]['weight']*self.length_multiplier(reaction)
        #limit based on maximum_reaction_number:        
        rankings = {k:rankings[k] for k in list(rankings)[0:maximum_reaction_number]}
        return(rankings)    

    def _get_weighted_reaction_rankings_2(
            self,
            weighted_random_compounds:dict,
            maximum_reaction_number:int = 20,
            shortest_path_method:str = 'Djikstra',
            ):
        """
        given a dictionary of weighted random compounds from self.get_weighted_random_compounds find the shortest path between: 
        1. compound 0 and compound 1 -> list
        2. add further reactions based on C0 and CN and C1 and CN 
        3. weight based on availability of how many weighted_random_compounds are present in the reaction
        4. weight based on edge weight and (optional) length multiplier for unreasonable large reactions (if option selected)
        returns a dictionary of ranked reactions and their weighting. 

        #TODO: check based on whether each element is present in the reaction. 
        """
        
        if len(weighted_random_compounds) == 1:
            return(None)

        c1 = weighted_random_compounds[0]
        c2 = weighted_random_compounds[1]
        # 1st find possible paths between compounds 1 & 2 and others
        possibilities = []
        for cn in weighted_random_compounds:
            try:
                possibilities.append(
                    [
                        x[1] for x in list(
                            nx.shortest_paths.all_shortest_paths(
                                G=self.graph, source=c1, target=cn)
                        )
                    ]
                )
                possibilities.append(
                    [
                        x[1] for x in list(
                            nx.shortest_paths.all_shortest_paths(
                                G=self.graph, source=c2, target=cn)
                        )
                    ]
                )
            except IndexError:
                pass

        possibilities = it.chain(*possibilities)
        #now weight by how many species in the reaction are in the weighted_random_compounds
        #more present compounds in a reaction is better 
        possibilities = {
            reaction_index:self.weight_filter(reaction_index,weighted_random_compounds) for reaction_index in possibilities
            }
        possibilities = dict(
            sorted(possibilities.items(),key=lambda item: item[1],reverse=True)#[0:maximum_reaction_number]
            )
        
        rankings = {}
        for i,reaction in enumerate(possibilities):
            try:
                weight = self.graph.get_edge_data(
                    u=c1,
                    v=reaction
                    )[0]['weight']*self.length_multiplier(reaction)
            except TypeError:
                weight = self.graph.get_edge_data(
                    u=c2,
                    v=reaction
                    )[0]['weight']*self.length_multiplier(reaction)
                    
                rankings[reaction] = weight 

        #limit based on maximum_reaction_number:
        #rankings = {k:v for i,(k,v) in enumerate(rankings.items()) if i <maximum_reaction_number}
        rankings = dict(sorted(rankings.items(),key=lambda item: item[1])[0:maximum_reaction_number])
        return(rankings)          

#done
###############################analysis.py

    def average_sampling(
        self,
            average: str = 'mean'
    ) -> dict:
        
        t = list(self.data)[0]
        p = list(self.data[t])[0]
        zeroth = list(self.data[t][p])[0]
        #if isinstance(xr[0],str):
        #    #data = str_to_int_dict(data)
        
        md = {}
        f_1 = {}
        
        for T in self.data:
            f_2 = {}
            m_2 = {}
            for P in self.data[T]:
                f_3 = {}
                m_3 = {}
                for c in self.data[T][P][zeroth]['data']:
                    f_4 = []
                    for x in self.data[T][P]:
                        if not x == zeroth:
                            f_4.append(self.data[T][P][x]['data'][c])
                    diff = [i - self.data[T][P][zeroth]['data'][c] for i in f_4]
                    f_3[c] = np.mean(f_4)/1e-6
                    m_3[c] = {'value':np.mean(diff)/1e-6,'variance':np.var(diff)/1e-6} # 2nd value is the variance and not the std deviation
                m_2[float(P)] = m_3
                f_2[float(P)] = f_3
            md[float(T)] = m_2
            f_1[float(T)] = f_2
        
        self.mean_data = md
        self.final_concs = f_1


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
    
    
