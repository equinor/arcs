

<img src="./static/ARCS_Logo.png" width="100" align="left" alt="ARCS Logo seasoned by ChatGPT"> 

<div id="toc">
  <ul style="list-style: none">
    <summary>
      <h2> <pre>ARCS - Automated Reactions for CO<sub>2</sub> Storage</pre> </h2>
    </summary>
  </ul>
</div>

Version 1.5.0

### Installation

```
git clone https://github.com/badw/arcs.git
cd arcs
pip install . 
```

### Examples

an example Jupyter Notebook can be found here: 


`./examples/example.ipynb`

an example DFT data is found here: 

`./app/data/dft_data.json`

which was run using VASP and the SCAN meta-GGA functional. 


simple usage: 

```
from arcs.generate import GraphGenerator
from arcs.traversal import Traversal
from arcs.generate import GenerateInitialConcentrations

graph = GraphGenerator().from_file(
    filename='../app/data/dft_data.json',
    temperature=248,
    pressure=20,
    max_reaction_length=5
)

concentrations = GenerateInitialConcentrations(graph=graph).update_ic(
    {'H2O':30,'O2':10,'SO2':10,'H2S':10,'NO2':10}
    )

t = Traversal(graph=graph)

results = t.sample(
  initial_concentrations=concentrations,
  ncpus=4,
  nsamples=1000
  )
```

`results` can then be analysed with `arcs.analysis.AnalyseSampling`

1. reaction statistics

```
from arcs.analysis import AnalyseSampling
import pandas as pd 

analysis = AnalyseSampling()
stats = pd.Series(analysis.reaction_statistics(data)).sort_values(ascending=False)
stats.head(10)
```

>```1 H2 + 1 SO2 = 1 O2 + 1 H2S              369
>1 H2O + 1 SO2 = 1 H2SO3                  270
>2 H2 + 1 O2 = 2 H2O                      227
>3 O2 + 2 H2S = 2 H2O + 2 SO2             163
>.
>.
>.
>```

2. Mean average data

```
average_data = pd.DataFrame(analysis.average_sampling(data))
average_data = average_data.loc[~(average_data==0).all(axis=1)]
average_data.sort_values(by='diff',inplace=True)
average_data.round(2)
```



>```
>compound  initial mean  diff  sem   std   var
>H2S        10.0   4.88 -5.12  0.10  4.75  22.53
>NO2        10.0   6.19 -3.81  0.10  4.85  23.48
>O2         10.0   6.24 -3.76  0.12  5.76  33.18
>.
>.
>.
>```

### ARCS App 

The app is created using the `plotly DASH` ([https://github.com/plotly/dash](https://github.com/plotly/dash)) framework.  

The `arcs-app` can be run from the `app` directory in the terminal through; 

```
cd app
python arcs-app.py
```


<p align="center">
 <img src="./static/ARCS-gui.png" height="300">
</p>
