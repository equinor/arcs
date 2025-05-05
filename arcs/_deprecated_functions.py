#class ApplyDataToReaction:
#    ''' this class applies the Gibbs data to a specific reaction'''
#    
#    def __init__(
#            self,
#            trange,
#            prange,
#            data,
#            nprocs
#            ):
#        self.trange = trange
#        self.prange = prange
#        reactions = data['reactions']
#        try:
#            self.reactions = {i:Equilibrium.from_string(r) for i,r in enumerate(reactions)}
#        except Exception:
#            self.reactions = {i:r for i,r in enumerate(reactions)}
#        self.compound_data = {k:data[k] for k in data.keys() if not k == 'reactions'}
#        self.nprocs = nprocs
#        self.barformat = '{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}'
#        
##    def _generate_data_serial(self,t,p): #serial
##        reactions = {i:{'e':r,
##            'k':ReactionGibbsandEquilibrium(t,p,self.compound_data).#equilibrium_constant(r),
##            'g':ReactionGibbsandEquilibrium(t,p,self.compound_data).#reaction_energy(r)} 
##                     for i,r in self.reactions.items()}
##        return(reactions)
#    
#    def _generate_data_serial(self,t,p):
#        rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
#        reactions = {i:rge.as_dict(r) for i,r in self.reactions.items()}
#        return(reactions)#

#    def generate_data(self,t,p): #multiprocessed#

#        manager = pmp.Manager()
#        queue = manager.Queue()
#        
#        def mp_function(reaction_keys,out_q):#

#            data = {}
#            for r in reaction_keys:
#                rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
#                data[r] = rge.as_dict(self.reactions[r])
#            out_q.put(data)#

#        resultdict = {}
#        r_keys = list(self.reactions.keys())
#        chunksize = int(math.ceil(len(self.reactions)/float(self.nprocs)))
#        processes = []#

#        for i in range(self.nprocs):
#            pr = pmp.Process(target=mp_function,
#                            args=(r_keys[chunksize*i:chunksize*(i+1)],queue))
#            processes.append(pr)
#            pr.start()#

#        for i in range(self.nprocs):
#            resultdict.update(queue.get(timeout=1800))#

#        for pr in processes:
#            pr.join()#

#        return(resultdict)#
#

#    def apply(self,serial=False):
#        data = {}
#        for t in self.trange:
#            pdat = {}
#            for p in self.prange:
#                if serial:
#                    pdat[p] = self._generate_data_serial(t,p)
#                else:
#                    pdat[p] = self.generate_data(t,p)
#            data[t] = pdat
#        self.data = data
#        return(self.data) 
#    
#    def save(self,filename='applied_reactions.p'):
#        pickle.dump(self.data,open(filename,'wb'))
#        print('data saved to: {}'.format(filename))