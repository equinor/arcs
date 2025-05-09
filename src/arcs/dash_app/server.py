from multiprocessing import Condition

from traitlets import default
from arcs.dash_app.domino import terminate_when_parent_process_dies
import dash
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
import dash_loading_spinners as dls
from dash import html
from dash import dash_table
from dash import dcc
from dash import ctx
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
from monty.serialization import loadfn
from arcs.generate import GenerateInitialConcentrations
from arcs.analysis import AnalyseSampling
from arcs.traversal import Traversal
from arcs.generate import GraphGenerator
from dash.exceptions import PreventUpdate


from arcs.dash_app.styling import keys_by_depth,_markdown_compound,format_concentrations_table_for_dash,format_statistics_table_for_dash,format_bar_chart

def start_dash(host: str, 
               port: int, 
               server_is_started: Condition,
): 
              # file_location):

    terminate_when_parent_process_dies()

    external_stylesheets = [dbc.themes.MINTY]
    
    load_figure_template("MINTY")

    
    graph = None 
    table4 = None
    #graph = html.P('None')
    #table4 = html.P('None')
    
    
    dft_filename= './data/dft_data.json'
    # generate graph from dft data

    #ambient_conditions = {'temperature':298, 'pressure':1}

    #default_concentrations = {'H2O':10,'SO2':10,'NO2':10}

    default_settings = {
        "exclude_co2": True,
        "max_compounds": 5,
        "discovery_threshold": 5,
        "maximum_reaction_number": 10,
        "max_steps": 5,
        "nsamples": 100,
        "ncpus": 4,
        "ceiling": 2000,
        "scale_largest": 0.2,
        "rank_small_reactions_higher":True,
        "rank_by_number_of_atoms":True,
        "shortest_path_method":'Djikstra'
    }


    #_graph = GraphGenerator().from_dict(
    #    dft_dict=default_dft_dict,
    #    temperature=ambient_settings['T'],
    #    pressure=ambient_settings['P'])

    ### generate Initial Concentrations
    #concs = GenerateInitialConcentrations(graph=_graph).update_ic(
    #    default_concentrations
    #    ,include_co2=False)



    
    ###################### layout of DASH template########################
    app = dash.Dash(
        __name__, 
        external_stylesheets=external_stylesheets
        )
    
    loading_spinner = dls.Rotate(
        id="loading-1",
        width=150,
        margin=9,
        speed_multiplier=0.8,
        color="rgba(58, 136, 254,1)",
        fullscreen=False,
        children=html.Div(id="loading-output-1"),
        fullscreen_style={"background-color": "rgba(0.1,0.1,0.1,0.2)"},
    )
    
    temperature_input = html.Div(
        [
            dbc.Label("Temperature (K)"),
            dbc.Input(
                placeholder="298",
                value=298,
                type="number",
                className="mb-3",
                id="temperature_input",
                debounce = True
                )
        ]
    )
    
    pressure_input = html.Div(
        [
            dbc.Label("Pressure (bar)"),
            dbc.Input(
                placeholder="1",
                value=1,
                type='number',
                className="mb-3",
                id='pressure_input',
                debounce = True
            )
        ]
    )
    
    concentrations_input = dbc.Stack(
        style={
            'textAlign': 'justified',
            "margin-left":"20px",
            "margin-right":"20px"
        },
        gap=3,
        children=[
            dash_table.DataTable(
                id='concentrations_input',
                columns=[{
                    'name': 'compound',
                    'id': 'index',
                    'editable': True
                },
                    {
                    'name': 'initial conc. (ppm)',
                    'id': 'initial',
                    'editable': True
                },
                ],
                data=[
                    {'index': 'H2O', 'initial': 30},
                    {'index': 'O2', 'initial': 10},
                    {'index': 'SO2', 'initial': 10},
                    {'index': 'NO2', 'initial': 0},
                    {'index':'H2S','initial':10}
                ],
                row_deletable=True,
                style_as_list_view=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    'padding-right': '30px',
                    'padding-left': '30px',
                    'text-align': 'center',
                    'marginLeft': 'auto',
                    'marginRight': 'auto'
                },
                style_table={
                    "overflow": "scroll",
                },
                fixed_rows={"headers": True},
    ),
    dbc.Button('add compound', id='addrows', n_clicks=0)
        ]
    )
    
    arcs_settings = dbc.Accordion(
        start_collapsed=True,
        children=[
            dbc.AccordionItem(
                title="Number of Samples",
                className="accordion",
                children=[
                    html.P(["Number of sampling events used to get a stochastic mean average.",html.Br(),"Default = 1000",html.Br(),"Recommended Amount > 500",html.Br()]),
                    dcc.Input(
                        id="nsamples",
                        value="1000",
                        debounce=True,
                        className="form-label mt-4",
                    ),
                ],
            ),
            dbc.AccordionItem(
                title="Maximum Random Walk Steps",
                className="accordion",
                children=[
                    html.P(["The maximum number of steps is the upper limit for how far down a reaction network a single sampling event goes.",html.Br(),"Default = 5",html.Br()]),
                    dcc.Input(
                        id="max_steps",
                        value="5",
                        debounce=True,
                        className="form-label mt-4",
                    ),
                ],
            ),
            dbc.AccordionItem(
                title="Discovery % Threshold",
                className="accordion",
                children=[
                    html.P(["The discovery % threshold, determines the amount (relative to the total concentration of compounds in the system) at which a certain compound is deemed 'chooseable' for ranking possible reactions. This is in order to weight the reactions in terms of the larger concentrations dominating the reaction probabilities.",html.Br(),"Default = 5%",html.Br()]),
                    dcc.Input(
                        id="discovery_threshold",
                        value="5",
                        debounce=True,
                        className="form-label mt-4",
                    ),
                ],
    
            ),
            dbc.AccordionItem(
                title="Concentration % Ceiling",
                className="accordion",
                children=[
                    html.P(["The concentration % ceiling determines the border by which a component should be scaled down to allow for reactions with smaller components. The component value is then scaled by the amount specified in 'Largest Concentrations Scale Value'. It is important to note that this is only for choosing reactions and not the final value which remains as the large amount. This helps to alleviate situations where no reactions are expected in ARCS due to overweighting towards the large concentration components.",html.Br(),"Default = 2000%",html.Br()])
                    ,
                    dcc.Input(
                        id="ceiling",
                        value="2000",
                        debounce=True,
                        className="form-label mt-4",
                    ),
                ],
            ),
            dbc.AccordionItem(
                title="Largest Concentrations Scale percentage",
                className="accordion",
                children=[
                    html.P(["Components which have reached the % ceiling in 'Concentration % Ceiling' are scaled by this amount.",html.Br(),"Default = 10% of original value",html.Br()]),
                    dcc.Input(
                        id="scale_largest",
                        value="10",
                        debounce=True,
                        className="form-label mt-4",
                    )
                ],
            ),
            dbc.AccordionItem(
                title="Max. Number of Reactions Considered in Choice",
                className="accordion",
                children=[
                    html.P(["The maximum number of equations suitable for a random choice after ranking all possible solutions.",html.Br(),"Default = 10",html.Br()]),
                    dcc.Input(
                        id="maximum_reaction_number",
                        value="10",
                        debounce=True,
                        className="form-label mt-4",
                    )
                ]
            ),
            dbc.AccordionItem(
                title="Max. Number of Compounds Considered in Choice",
                className="accordion",
                children=[
                    html.P(["This corresponds to the largest number of compounds that can be considered for ranking at any one time.",html.Br(),"Default = 5",html.Br()]),  
                    dcc.Input(
                        id="max_compounds",
                        value="5",
                        debounce=True,
                        className="form-label mt-4",
                    )
                ]
            ),
            dbc.AccordionItem(
                title="Shortest Path Method",
                className="accordion",
                children=[
                    html.P(["There are multiple algorithms that can be used to find the shortest path between two components in the reaction graph.",html.Br(),"Implemented in this work are the 'Bellman-Ford' and 'Dijkstra' methods.",html.Br(),"Default = Dijkstra",html.Br()]),
                    dbc.RadioItems(
                        id="shortest_path_method",
                        className="btn btn-outline-primary",
                        options=[
                            {"label": "Bellman-Ford", "value": "Bellman-Ford"},
                            {"label": "Dijkstra", "value": "Dijkstra"},
                        ],
                        value="Dijkstra",
                    )
                ]
            ),

            dbc.AccordionItem(
                title="Rank Smaller Reactions Higher",
                className="accordion",
                children=[
                    html.P(
                        ["Sometimes the simplest solutions are the most plausible. Rank Smaller reactions higher in the list.",html.Br(), "Default = True",html.Br()]),
                    dbc.RadioItems(
                        id="rank_small_reactions_higher",
                        className="btn btn-outline-primary",
                        options=[
                            {"label":"True","value":True},
                            {"label":"False","value":False},
                        ],
                        value=True,
                    )
                ]
            ),
        ],
    ),
    
    submit_button = dbc.Button(
                children="Run",
                id="submit-val",
                n_clicks=0,
                className="btn btn-success",
                style={'float': 'left',"margin-right":"1rem"}
            )
    
    offcanvas = html.Div(
        style={
            'textAlign': 'justified',
            "margin-left": "1px",
            "margin-right": "1px",
        },
        children=[
            dbc.Button("Settings", id="open-offcanvas", n_clicks=0,className='btn btn-info',style={'float': 'left',"margin-right":"1rem"}),
            dbc.Offcanvas(
                children=[
                    dbc.Stack(
                        [
                    dbc.Card(
                        [
                            dbc.CardHeader("ARCS Settings"),
                            dbc.CardBody(arcs_settings)
                        ],
                        #color='dark',
                    ),
                        ],
                        gap=3
                    )
                ],
                id="offcanvas",
                is_open=False,
                scrollable=True,
                style={"width":"50rem"}
            )
        ]
    )
    
    most_frequent_reactions = html.Div(
        id="reaction-stats",
        style={"align": "center"},
        children=table4,
    )
    
    logos = html.Div( # needs work
        children=[
            html.Img(
                src=app.get_asset_url(
                    "images/logos.png"
                ),
                style={
                    "width": "10%",
                    "height": "10%",
                    "padding": "0.05rem",
                    "align": "end",
                },
            ),
        ],
    ),
    
    results_concentration_table = html.Div(
        id="final_concs_table",
        children=None,
    ),

    results_concentration_bar_chart = html.Div(
        id="final_concs_barchart",
        children=graph
    ),
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################### layout
    
    app.layout = html.Div(
        style={'padding': '5rem'},
        children=[
            dbc.Row(dbc.Col(logos)),
            loading_spinner,
            dbc.Row(
                [
                    html.P("ARCS 1.5.0"),
                    html.H3(["Automated Reactions for ","C", "O", html.Sub(2), " Conversion (ARCS)"]),
                    html.Div(
                        [offcanvas,
                         submit_button,
                         ]
                        ),
                    html.Div(
                        # for updating the concentrations to be used in ARCS (no need for displaying)
                        id="placeholder1",
                        children=None,
                        style={"display": "none"},
                    ),
                    html.Div(
                        id="placeholder2",  # for updating the settings to be used in ARCS
                        children=None,
                        style={"display": "none"},
                    ),
                    html.Div(
                        id="placeholder3",  # placeholder3 is for updating the sliders data to be used in ARCS
                        children=None,
                        style={"display": "none"},
                    ),
                    html.Div(
                        id="placeholder4",  # placeholder4 is for updating the upload file 
                        children=None,
                        style={"display": "none"},
                    ),
                ]
            ),
            dbc.Tabs(
                style={'padding':'2rem','align':'center'},
                children=[
                    dbc.Tab(
                        className='nav nav-tabs',
                        label='Inputs',
                        children=[
                            dbc.Col(
                                children=[dbc.Stack(
                                    gap=3,
                                    children=[
                                        dbc.Card(
                                            children=[
                                                dbc.CardHeader('Input Concentrations'),
                                                dbc.CardBody(concentrations_input),
                                            ],
                                            #color='dark'
                                        ),
                                        dbc.Card(
                                            children=[
                                                dbc.CardHeader('Conditions'),
                                                dbc.CardBody([temperature_input,
                                                              pressure_input])
                                                #dbc.CardBody(sliders),
                                            ],
                                            #color='dark'
                                        )
                                    ]
                                )
                                ]
                            )
                        ]
                    ),
                    dbc.Tab(
                        label='Output Concentrations',
                        children=[
                            dbc.Stack(
                                gap=3,
                                children=[
                                    dbc.Card(
                                        children=[
                                            dbc.CardHeader(
                                                "Change in Concentrations"),
                                            dbc.CardFooter(
                                                dbc.Tabs(
                                                    style={'padding':'2rem'},
                                                    children=[
                                                        dbc.Tab(
                                                            results_concentration_bar_chart, label='BarChart'),
                                                        dbc.Tab(
                                                            results_concentration_table, label='Table')
                                                    ],
                                                )
                                            ),
                                        ],
                                        #color='dark'
                                    ),
                                ]
                            )
                        ]
                    ),
                    dbc.Tab(
                        label='Reactions',
                        children=[
                            dbc.Stack(
                                gap=3,
                                children=[
                                    dbc.Card(
                                        children=[
                                            dbc.CardHeader(
                                                'Most Frequent Reactions'),
                                            dbc.CardBody(
                                                most_frequent_reactions),
                                        ],
                                        #color='dark'
                                    ),
                                ]
                            )
                        ]
                        ),
                        ]
    
                        ),
                        ]
                    )
    
    
    
    
    
    
    
    
    
    
    
    
    ###############################################################################################    
    #################app callbacks
    #off canvas
    @app.callback(
        Output("offcanvas", "is_open"),
        Input("open-offcanvas", "n_clicks"),
        [State("offcanvas", "is_open")],
    )
    def toggle_offcanvas(n1, is_open):
        if n1:
            return not is_open
        return is_open
    
    #update concentrations table (new!)
    @app.callback(
            Output('concentrations_input','data'),
            Input('addrows','n_clicks'),
            State('concentrations_input','data'),
            State('concentrations_input','columns'))
    def add_row(n_clicks,rows,columns):
        if n_clicks > 0:
            rows.append({c['id']: '' for c in columns})
        else:
            raise PreventUpdate
        return(rows)
    
    ####update T and P and load a new graph
    @app.callback(
        Output("placeholder3", "children"),
        [
            Input("temperature_input", "value"),
            Input("pressure_input", "value"),
        ],
    )
    def update_t_and_p(*inputs):
        global ambient_conditions
        ambient_conditions = {
            'temperature':float(inputs[0]),
            'pressure':float(inputs[1])
            }
    
    ####update the concentrations
    @app.callback(
            Output('placeholder1','children'),
            Input('concentrations_input','data'),
            Input('concentrations_input','columns'),
    )
    def update_concentrations(rows,columns):
        global default_concentrations
        default_concentrations = {
                row.get('index',None):float(row.get('initial',None)) for row in rows
                }

    # update settings
    @app.callback(
        Output("placeholder2", "children"),
        [
            Input("nsamples", "value"),
            Input("max_steps", "value"),
            Input("discovery_cutoff", "value"),
            Input("ceiling", "value"),
            Input("scale_largest", "value"),
            Input("maximum_reaction_number", "value"),
            Input("max_compounds", "value"),
            Input("shortest_path_method", "value"),
            Input("rank_small_reactions_higher","value")
        ],
    )
    def update_settings(inputs):
        default_settings["nsamples"]=int(inputs[1])
        default_settings["max_steps"]=int(inputs[2])
        default_settings["discovery_threshold"]=float(inputs[3]) / 100
        default_settings["ceiling"]=int(inputs[4])
        default_settings["scale_largest"]=float(inputs[5])
        default_settings["maximum_reaction_number"]=int(inputs[6])
        default_settings["max_compounds"]=int(inputs[7])
        default_settings["shortest_path_method"]=str(inputs[8])
        default_settings["rank_small_reactions_higher"]=bool(inputs[10])

    @app.callback(
        [
            Output("reaction-stats", "children"),
            Output("final_concs_table", "children"), 
            Output("final_concs_barchart", "children"),
            Output("loading-output-1", "children"),
        ],
        Input("submit-val", "n_clicks"),
    )
    def apprun(btn1):
        global default_concentrations
        global ambient_conditions

        if "submit-val" == ctx.triggered_id:
            graph = GraphGenerator().from_file(
                filename=dft_filename,
                temperature=ambient_conditions['temperature'],
                pressure=ambient_conditions['pressure'],
                max_reaction_length=3 # for quick testing and debugging - can make this a setting later. 
            )

            default_concentrations = GenerateInitialConcentrations(graph=graph).update_ic(update_dict=default_concentrations,include_co2=True)
            
            t = Traversal(graph=graph)
            
            print(ambient_conditions)
            print(default_concentrations)
            results = t.sample(
                initial_concentrations = default_concentrations,
                ncpus=default_settings['ncpus'],
                nsamples=default_settings['nsamples'],
                tqdm_kws={'disable':False} #can turn off for debugging
            )
            
            analysis = AnalyseSampling()

            reaction_statistics = pd.Series(
                analysis.reaction_statistics(results)
                ).sort_values(ascending=False)
            
            reaction_statistics = pd.DataFrame(
                reaction_statistics,columns=['frequency']
                ).reset_index()
            
            formatted_stats_table = format_statistics_table_for_dash(reaction_statistics)
            
            average_data = pd.DataFrame(analysis.average_sampling(results))
            average_data = average_data.loc[~(average_data==0).all(axis=1)]
            average_data.sort_values(by='diff',inplace=True)
            

            formatted_average_table = format_concentrations_table_for_dash(average_data.round(2))

            formatted_bar_chart = format_bar_chart(average_data)

            return (
                [formatted_stats_table],
                [formatted_average_table],
                [formatted_bar_chart],
                [None]
            )



    with server_is_started:
        server_is_started.notify()
    
    app.run(debug=False, port=port, host=host)





#    
#            return [
#                [metadata_table],
#                [stats_table],
#                [paths_table],
#                [diff_table],
#                [resultsgraph],
#                [None],
#            ]
    

