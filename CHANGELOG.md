ARCS 1.5.0 

`setup_functions.py`

- removed ApplyDataToReactions as this class was now defunct 

`ReactionGibbsandEquilibrium`
- removed a lot of chempy dependency to make it faster 
- implemented reactit output which is more pythonic and faster 
- docstrings added 
- code readability now much better 
- tests added for this class

`GraphGenerator` 
- removed trange and prange as it is not necessary anymore 
- removed multiprocessing functions as they were not necessary anymore 
- much faster
- code readability much better
- added tests for this class