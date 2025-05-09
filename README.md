

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
>3 H2 + 1 SO2 = 2 H2O + 1 H2S             136
>1 H2O + 1 NO2 + 1 NO = 2 HNO2            115
>1 H2 + 1 SO2 + 1 NO2 = 1 H2SO3 + 1 NO    114
>1 H2O + 2 NO2 = 1 HNO3 + 1 HNO2           78
>1 H2 + 1 H2SO4 = 2 H2O + 1 SO2            72
>1 H2 + 1 NO2 = 1 H2O + 1 NO               68
>```

2. Mean average data

```
average_data = pd.DataFrame(analysis.average_sampling(data))
average_data = average_data.loc[~(average_data==0).all(axis=1)]
average_data.sort_values(by='diff',inplace=True)
average_data.round(2)
```

>```        initial   mean  diff   sem   std    var
>H2S        10.0   4.88 -5.12  0.10  4.75  22.53
>NO2        10.0   6.19 -3.81  0.10  4.85  23.48
>O2         10.0   6.24 -3.76  0.12  5.76  33.18
>S8          0.0   0.07  0.07  0.01  0.26   0.07
>NH3         0.0   0.12  0.12  0.02  0.85   0.73
>NOHSO4      0.0   0.17  0.17  0.02  0.94   0.89
>HNO3        0.0   0.45  0.45  0.03  1.52   2.30
>HNO2        0.0   0.48  0.48  0.04  1.78   3.17
>H2O        30.0  30.48  0.48  0.11  5.51  30.31
>H2SO3       0.0   0.51  0.51  0.04  2.17   4.71
>N2          0.0   0.54  0.54  0.03  1.53   2.33
>H2SO4       0.0   0.65  0.65  0.05  2.47   6.08
>NO          0.0   1.45  1.45  0.07  3.45  11.89
>H2          0.0   2.19  2.19  0.08  4.04  16.32
>SO2        10.0  12.65  2.65  0.11  5.55  30.82

### ARCS App 

The app is created using the `plotly DASH` ([https://github.com/plotly/dash](https://github.com/plotly/dash)) framework.  

The `arcs-app` can be run from the `app` directory in the terminal through; 

```
python app/arcs-app.py
```


<p align="center">
 <img src="./static/ARCS-gui.png" width="300" height="300">
</p>
