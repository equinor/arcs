from multiprocessing import Condition
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
from arcs.setup_functions import GenerateInitialConcentrations
from arcs.analysis import AnalyseSampling
from arcs.traversal import Traversal
import pickle
import warnings


def start_dash(host: str, port: int, server_is_started: Condition, file_location):
    terminate_when_parent_process_dies()
    external_stylesheets = [dbc.themes.QUARTZ]

    load_figure_template("QUARTZ")

    # data and sliders

    def load_data(filename):  # has to be a .json file in dict format
        data = loadfn(filename)
        return data

    def keys_by_depth(dict_, depth=0, output=None):
        if output is None:
            output = {}
        if depth not in output:
            output[depth] = set()
        for key in dict_:
            output[depth].add(key)
            if isinstance(dict_[key], dict):
                keys_by_depth(dict_[key], depth + 1, output)
        return output

    def _markdown_compound(_string):
        md = []
        for i in _string:
            try:
                int(i)
                md.append("<sub>{}</sub>".format(int(i)))
            except ValueError:
                md.append(i)
        return "".join(md)

    def make_sliders(import_data, labels):
        def sliderform(key, slider, label):
            # this is very janky - need to find a better way of storing the data
            correct_order = sorted([float(x) for x in list(slider)])
            # marks = {float(x):str('{:.2F}'.format(float(x))) for x in list(correct_order)}
            marks = {
                int(x): {"label": str(x), "style": {"color": "rgba(0.1,0.1,0.1,0)"}}
                for x in correct_order
            }
            minval = float(list(correct_order)[0])
            maxval = float(list(correct_order)[-1])

            return html.Div(
                style={"padding": "2rem"},
                children=[
                    # dbc.Col(children=html.Label(
                    #        children=label[key]), width=2),
                    dbc.Col(
                        children=[
                            html.Label(children=label[key]),
                            dcc.Slider(
                                className="form-range",
                                id="slider-{}".format(key),
                                min=minval,
                                max=maxval,
                                step=None,
                                marks=marks,
                                value=minval,
                                updatemode="drag",
                                tooltip={
                                    "placement": "bottom",
                                    "always_visible": True,
                                    "style": {
                                        "color": "LightSteelBlue",
                                        "fontSize": "20px",
                                    },
                                },
                            ),
                        ]
                        # width=10,
                    ),
                ],
            )

        slider_keys = keys_by_depth(import_data)
        sliders = []
        for k in slider_keys:
            sliders.append(sliderform(k, sorted(slider_keys[k]), labels))
        return sliders

    # run data fields
    g = pickle.load(open(file_location + "SCAN_graph.p", "rb"))
    t = Traversal(graph=g, reactions=file_location + "SCAN_reactions.p")

    graph = dbc.Alert("No Data", color="light")  # None #html.P('None')
    table4 = dbc.Alert("No Data", color="light")  # None #html.P('None')
    table5 = dbc.Alert("No Data", color="light")  # None #html.P('None')

    meta = dbc.Alert("Data Shown When Run", color="secondary")  # None #html.P('None')
    sliders = make_sliders(g, labels={0: "T (K)", 1: "P (bar)"})

    gic = GenerateInitialConcentrations(g)
    gic.all_zero(include_co2=False)
    gic.update_ic({"SO2": 10e-6, "NO2": 50e-6, "H2S": 30e-6, "H2O": 20e-6})
    concs = gic.ic
    settings = {
        "nprocs": 1,
        "sample_length": 320,
        "max_rank": 10,
        "max_compounds": 5,
        "probability_threshold": 0.1,
        "path_depth": 5,
        "ceiling": 2000,
        "scale_highest": 0.2,
        "rank_small_reactions_higher": True,
    }
    ambient_settings = {"T": None, "P": None}

    ###################### layout of DASH template########################
    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    loading_spinner = dls.Rotate(
        id="loading-1",
        width=150,
        margin=9,
        speed_multiplier=0.8,
        color="rgba(58, 136, 254,1)",
        fullscreen=True,
        children=html.Div(id="loading-output-1"),
        fullscreen_style={"background-color": "rgba(0.1,0.1,0.1,0.2)"},
    )

    concentrations_table = dbc.Stack(
        style={"textAlign": "justified", "margin-left": "20px", "margin-right": "20px"},
        gap=3,
        children=[
            dash_table.DataTable(
                id="concentrations_table",
                columns=[
                    {"name": "compound", "id": "index", "editable": True},
                    {"name": "initial conc. (ppm)", "id": "initial", "editable": True},
                ],
                data=[
                    {"index": "H2O", "initial": 30},
                    {"index": "O2", "initial": 10},
                    {"index": "SO2", "initial": 10},
                    {"index": "NO2", "initial": 0},
                    {"index": "H2S", "initial": 10},
                ],
                row_deletable=True,
                style_as_list_view=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "30px",
                    "padding-left": "30px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_table={
                    "overflow": "scroll",
                },
                fixed_rows={"headers": True},
            ),
            dbc.Button("add compound", id="addrows", n_clicks=0),
        ],
    )

    arcs_settings = (
        dbc.Accordion(
            start_collapsed=True,
            children=[
                dbc.AccordionItem(
                    title="Samples",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "Number of sampling events used to get a stochastic mean average.",
                                html.Br(),
                                "Default = 100",
                                html.Br(),
                                "Recommended Amount > 500",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="samples",
                            value="100",
                            debounce=True,
                            # style={
                            # "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Maximum Path Depth",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "The maximum path depth is the upper limit for how far down a reaction network a single sampling event goes.",
                                html.Br(),
                                "Default = 5",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="pathdepth",
                            value="5",
                            debounce=True,
                            # style={
                            #    "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Discovery % Cutoff",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "The discovery % cutoff, determines the amount (relative to the total concentration of compounds in the system) at which a certain compound is deemed 'chooseable' for ranking possible reactions. This is in order to weight the reactions in terms of the larger concentrations dominating the reaction probabilities.",
                                html.Br(),
                                "Default = 5%",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="probability_cutoff",
                            value="5",
                            debounce=True,
                            # style={
                            #    "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Concentration % Ceiling",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "The concentration % ceiling determines the border by which a component should be scaled down to allow for reactions with smaller components. The component value is then scaled by the amount specified in 'Largest Concentrations Scale Value'. It is important to note that this is only for choosing reactions and not the final value which remains as the large amount. This helps to alleviate situations where no reactions are expected in ARCS due to overweighting towards the large concentration components.",
                                html.Br(),
                                "Default = 500%",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="ceiling",
                            value="500",
                            debounce=True,
                            # style={
                            #    "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Largest Concentrations Scale Value",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "Components which have reached the % ceiling in 'Concentration % Ceiling' are scaled by this amount.",
                                html.Br(),
                                "Default = 0.1",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="scale_highest",
                            value="0.1",
                            debounce=True,
                            # style={
                            #    "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Max. Number of Reactions Considered in Choice",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "The maximum number of equations suitable for a random choice after ranking all possible solutions.",
                                html.Br(),
                                "Default = 5",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="max_rank",
                            value="5",
                            debounce=True,
                            # style={
                            #    "backgroundColor": backgroundcolours,
                            #    "color": "white",
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Max. Number of Compounds Considered in Choice",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "This corresponds to the largest number of compounds that can be considered for ranking at any one time.",
                                html.Br(),
                                "Default = 5",
                                html.Br(),
                            ]
                        ),
                        dcc.Input(
                            id="max_compounds",
                            value="5",
                            debounce=True,
                            # style={
                            #    "color": "white",
                            #    "backgroundColor": backgroundcolours,
                            # },
                            className="form-label mt-4",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Graph Sampling Method",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "There are multiple algorithms that can be used to find the shortest path between two components in the reaction graph.",
                                html.Br(),
                                "Implemented in this work are the 'Bellman-Ford' and 'Dijkstra' methods.",
                                html.Br(),
                                "Default = Dijkstra",
                                html.Br(),
                            ]
                        ),
                        dbc.RadioItems(
                            id="method",
                            className="btn btn-outline-primary",
                            options=[
                                {"label": "Bellman-Ford", "value": "Bellman-Ford"},
                                {"label": "Dijkstra", "value": "Dijkstra"},
                            ],
                            value="Dijkstra",
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Include CO2 as a reactant?",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "Including CO",
                                html.Sub(2),
                                " as a reactant or assume that CO",
                                html.Sub(2),
                                " is a background solvent.",
                                html.Br(),
                                "Default = False",
                                html.Br(),
                            ]
                        ),
                        dbc.RadioItems(
                            id="include_co2",
                            className="btn btn-outline-primary",
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=False,
                        ),
                    ],
                ),
                dbc.AccordionItem(
                    title="Rank Smaller Reactions Higher",
                    className="accordion",
                    children=[
                        html.P(
                            [
                                "Sometimes the simplest solutions are the most plausible. Rank Smaller reactions higher in the list.",
                                html.Br(),
                                "Default = True",
                                html.Br(),
                            ]
                        ),
                        dbc.RadioItems(
                            id="rank_small_reactions",
                            className="btn btn-outline-primary",
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                        ),
                    ],
                ),
            ],
        ),
    )

    submit_button = dbc.Button(
        children="Run",
        id="submit-val",
        n_clicks=0,
        className="btn btn-success",
        style={"float": "left", "margin-right": "1rem"},
    )

    metadatatable = html.Div(
        id="metadata",
        style={"align": "center"},
        children=meta,
        # className="table table-secondary",
    )

    offcanvas = html.Div(
        style={
            "textAlign": "justified",
            "margin-left": "1px",
            "margin-right": "1px",
            # "padding":"2px",
        },
        children=[
            dbc.Button(
                "Settings",
                id="open-offcanvas",
                n_clicks=0,
                className="btn btn-info",
                style={"float": "left", "margin-right": "1rem"},
            ),
            dbc.Offcanvas(
                children=[
                    dbc.Stack(
                        [
                            dbc.Card(
                                [
                                    dbc.CardBody(arcs_settings),
                                    dbc.CardFooter("ARCS Settings"),
                                ]
                            ),
                            dbc.Card(
                                [
                                    dbc.CardBody(metadatatable),
                                    dbc.CardFooter("System Data"),
                                ]
                            ),
                        ],
                        gap=3,
                    )
                ],
                id="offcanvas",
                is_open=False,
                scrollable=True,
                style={"width": "50rem"},
            ),
        ],
    )

    most_frequent_reactions = html.Div(
        id="reaction-stats",
        style={"align": "center"},
        children=table4,
        # className="table table-primary",
    )

    most_frequent_paths = (
        html.Div(
            id="reaction-paths",
            style={"align": "center"},
            children=table5,
            # className="table table-primary",
        ),
    )

    logos = (
        html.Div(  # needs work
            children=[
                html.Img(
                    src=app.get_asset_url("images/logos.png"),
                    style={
                        # "width": "50%",
                        "height": "30%",
                        "padding": "0.05rem",
                        "align": "end",
                    },
                ),
            ],
        ),
    )

    results_concentration_graph = (
        html.Div(
            id="final_concs",
            children=None,
        ),
    )

    finalgraph = (html.Div(id="output-graph", children=graph),)

    ################################### layout

    app.layout = html.Div(
        style={"padding": "5rem"},
        children=[
            dbc.Row(dbc.Col(logos)),
            dbc.Row(
                [
                    html.P("ARCS 1.4.0"),
                    html.H3(
                        [
                            "Automated Reactions for ",
                            "C",
                            "O",
                            html.Sub(2),
                            " Conversion (ARCS)",
                        ]
                    ),
                    html.Div([offcanvas, submit_button]),
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
                ]
            ),
            loading_spinner,
            dbc.Tabs(
                style={"padding": "2rem", "align": "center"},
                children=[
                    dbc.Tab(
                        className="nav nav-tabs",
                        label="Inputs",
                        children=[
                            dbc.Col(
                                children=[
                                    dbc.Stack(
                                        gap=3,
                                        children=[
                                            dbc.Card(
                                                children=[
                                                    dbc.CardHeader(
                                                        "Input Concentrations"
                                                    ),
                                                    dbc.CardBody(concentrations_table),
                                                ]
                                            ),
                                            dbc.Card(
                                                children=[
                                                    dbc.CardHeader("Conditions"),
                                                    dbc.CardBody(sliders),
                                                ]
                                            ),
                                        ],
                                    )
                                ]
                            )
                        ],
                    ),
                    dbc.Tab(
                        label="Output Concentrations",
                        children=[
                            dbc.Stack(
                                gap=3,
                                children=[
                                    dbc.Card(
                                        children=[
                                            dbc.CardHeader("Change in Concentrations"),
                                            dbc.CardFooter(
                                                dbc.Tabs(
                                                    style={"padding": "2rem"},
                                                    children=[
                                                        dbc.Tab(
                                                            finalgraph, label="Graph"
                                                        ),
                                                        dbc.Tab(
                                                            results_concentration_graph,
                                                            label="Table",
                                                        ),
                                                    ],
                                                )
                                            ),
                                        ]
                                    ),
                                ],
                            )
                        ],
                    ),
                    dbc.Tab(
                        label="Reactions",
                        children=[
                            dbc.Stack(
                                gap=3,
                                children=[
                                    dbc.Card(
                                        children=[
                                            dbc.CardBody(most_frequent_reactions),
                                            dbc.CardHeader("Most Frequent Reactions"),
                                        ]
                                    ),
                                    dbc.Card(
                                        children=[
                                            dbc.CardBody(most_frequent_paths),
                                            dbc.CardHeader(
                                                "Most Frequent Paths (alpha)"
                                            ),
                                        ]
                                    ),
                                ],
                            )
                        ],
                    ),
                ],
            ),
        ],
    )

    #####################app callbacks
    # off canvas
    @app.callback(
        Output("offcanvas", "is_open"),
        Input("open-offcanvas", "n_clicks"),
        [State("offcanvas", "is_open")],
    )
    def toggle_offcanvas(n1, is_open):
        if n1:
            return not is_open
        return is_open

    # update concentrations table (new!)
    @app.callback(
        Output("concentrations_table", "data"),
        Input("addrows", "n_clicks"),
        State("concentrations_table", "data"),
        State("concentrations_table", "columns"),
    )
    def add_row(n_clicks, rows, columns):
        if n_clicks > 0:
            rows.append({c["id"]: "" for c in columns})
        return rows

    # update the concentrations
    @app.callback(
        Output("placeholder1", "children"),
        Input("concentrations_table", "data"),
        Input("concentrations_table", "columns"),
    )
    def update_concentrations(rows, columns):
        for k, v in concs.items():  # reset values
            if not k == "CO2":
                concs[k] = 0
        for row in rows:
            spec = row.get("index", None)
            num = row.get("initial", None)
            if spec in list(concs.keys()):
                concs[spec] = float(num) * 1e-6
            else:
                pass

    # update settings
    @app.callback(
        Output("placeholder2", "children"),
        [
            Input("samples", "value"),
            Input("pathdepth", "value"),
            Input("probability_cutoff", "value"),
            Input("ceiling", "value"),
            Input("scale_highest", "value"),
            Input("max_rank", "value"),
            Input("max_compounds", "value"),
            Input("method", "value"),
            Input("include_co2", "value"),
            Input("rank_small_reactions", "value"),
        ],
    )
    def update_settings(*inputs):
        settings["sample_length"] = int(inputs[0])
        settings["path_depth"] = int(inputs[1])
        settings["probability_threshold"] = float(inputs[2]) / 100
        settings["ceiling"] = int(inputs[3])
        settings["scale_highest"] = float(inputs[4])
        settings["max_rank"] = int(inputs[5])
        settings["max_compounds"] = int(inputs[6])
        settings["method"] = str(inputs[7])
        settings["include_co2"] = bool(inputs[8])
        settings["rank_small_reactions_higher"] = bool(inputs[9])

    # update T and P
    @app.callback(
        Output("placeholder3", "children"),
        [Input("slider-0", "value"), Input("slider-1", "value")],
    )
    def update_t_and_p(*inputs):
        ambient_settings["T"] = int(inputs[0])
        ambient_settings["P"] = int(inputs[1])

    @app.callback(
        [
            Output("metadata", "children"),
            Output("reaction-stats", "children"),
            Output("reaction-paths", "children"),
            Output("final_concs", "children"),
            Output("output-graph", "children"),
            Output("loading-output-1", "children"),
        ],
        Input("submit-val", "n_clicks"),
    )
    def apprun(btn1):
        if "submit-val" == ctx.triggered_id:
            warnings.simplefilter("ignore")
            t.run(
                trange=[ambient_settings["T"]],
                prange=[ambient_settings["P"]],
                save=False,
                ic=concs,
                **settings,
            )

            metadata = t.metadata
            metadata = pd.Series(metadata).reset_index()
            metadata = metadata.rename(columns={0: "value"})

            metadata_table = dash_table.DataTable(
                columns=[
                    {"name": i, "id": i, "type": "text", "presentation": "markdown"}
                    for i in metadata.columns
                ],
                data=metadata.to_dict("records"),
                style_as_list_view=False,
                cell_selectable=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "30px",
                    "padding-left": "30px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_table={
                    "overflow": "scroll",
                },
                markdown_options={"html": True, "link_target": "_self"},
                # fixed_rows={"headers": True},
            )
            #####updating concentrations table
            df_d = (
                pd.DataFrame(
                    t.initfinaldiff[ambient_settings["T"]][ambient_settings["P"]]
                )
                .round(1)
                .drop("CO2")
            )
            df_d = df_d.T
            df_d = df_d.T.rename({x: _markdown_compound(x) for x in df_d})
            df_d = df_d.reset_index()
            df_d = df_d.rename(
                {
                    "index": "compound",
                    "initial": "initial (ppm)",
                    "final": "final (ppm)",
                    "change": "change (ppm)",
                }
            )

            diff_table = dash_table.DataTable(
                columns=[
                    {"name": i, "id": i, "type": "text", "presentation": "markdown"}
                    for i in df_d.columns
                ],
                data=df_d.to_dict("records"),
                style_as_list_view=True,
                cell_selectable=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "30px",
                    "padding-left": "30px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_table={
                    "overflow": "scroll",
                },
                # fixed_rows={"headers": True},
                # style_cell_conditional=[
                #    {"if": {"column_id": "index"},
                #        "width": "10%", "textAlign": "left"},
                #    {
                #        "if": {"column_id": ["initial", "final", "change"]},
                #        "width": "10%",
                #        "textAlign": "left",
                #    },
                # ],
                markdown_options={"html": True, "link_target": "_self"},
            )

            ### statistics
            analyse = AnalyseSampling(t.data, markdown=True)
            analyse.reaction_statistics()
            analyse.mean_sampling()
            analyse.reaction_paths()

            mean = analyse.mean_data[ambient_settings["T"]][ambient_settings["P"]]

            paths = pd.DataFrame(
                analyse.common_paths[ambient_settings["T"]][ambient_settings["P"]]
            )
            stats = pd.DataFrame(
                analyse.stats[ambient_settings["T"]][ambient_settings["P"]]
            )
            stats_table = dash_table.DataTable(
                columns=[
                    {
                        "name": "Reactions",
                        "id": "index",
                        "type": "text",
                        "presentation": "markdown",
                    },
                    {
                        "name": "k",
                        "id": "k",
                        "type": "text",
                        "presentation": "markdown",
                    },
                    {
                        "name": "Frequency",
                        "id": "frequency",
                        "type": "text",
                        "presentation": "markdown",
                    },
                ],
                data=stats.to_dict("records"),
                style_as_list_view=False,
                cell_selectable=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "10px",
                    "padding-left": "10px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_header={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "10px",
                    "padding-left": "10px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_table={
                    "overflow": "scroll",
                },
                # fixed_rows={"headers": True},
                markdown_options={"html": True, "link_target": "_self"},
            )

            paths_table = dash_table.DataTable(
                columns=[
                    {
                        "name": "Paths",
                        "id": "paths",
                        "type": "text",
                        "presentation": "markdown",
                    },
                    {
                        "name": "k",
                        "id": "k",
                        "type": "text",
                        "presentation": "markdown",
                    },
                    {
                        "name": "Frequency",
                        "id": "frequency",
                        "type": "text",
                        "presentation": "markdown",
                    },
                ],
                # columns=[{'name':i,'id':i} for i in new_stats.columns],
                data=paths.to_dict("records"),
                style_as_list_view=False,
                cell_selectable=False,
                style_cell={
                    "font_family": "helvetica",
                    "align": "center",
                    "padding-right": "10px",
                    "padding-left": "10px",
                    "text-align": "center",
                    "marginLeft": "auto",
                    "marginRight": "auto",
                },
                style_table={
                    "overflow": "scroll",
                },
                # fixed_rows={"headers": True},
                markdown_options={"html": True, "link_target": "_self"},
            )
            df_m_t = pd.DataFrame(mean).T
            df_m_t = df_m_t[df_m_t["value"] != 0]

            df_m = pd.DataFrame(
                {
                    "comps": list(df_m_t.T.keys()),
                    "values": df_m_t["value"].values,
                    "variance": df_m_t["variance"].values,
                    "variance_minus": -df_m_t["variance"].values,
                }
            )
            # maxval=np.max(
            #    [np.abs(df_m["values"].min()), np.abs(df_m["values"].max())]
            # )
            # ymin, ymax=[-maxval, maxval]

            fig = px.bar(
                df_m,
                x="comps",
                y="values",
                error_y="variance",
                error_y_minus="variance_minus",
                labels={"comps": "", "values": "\u0394 ppm"},
                color="values",
                color_continuous_scale="tropic_r",
                hover_data={
                    "values": False,
                    "comps": False,
                    "variance": False,
                    "error": (":.2e", df_m["variance"]),
                    "specie": df_m["comps"],
                    "PPM": (":.1f", df_m["values"]),
                },
                # width=500,height=500
            )
            fig.update_layout(
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                hovermode="closest",
                hoverlabel=dict(font_size=16),
                coloraxis_showscale=False,
            )
            fig.update_xaxes(showgrid=False, tickangle=-60, tickmode="linear")
            # try:
            #    dtick=int(int(ymax - ymin) / 10)
            # except:
            #    dtick = None
            #    pass
            # fig.update_yaxes(
            #    showgrid=True,
            #    tickmode="linear",
            #    range=[ymin - 2, ymax + 2],
            #    dtick=dtick,
            # )

            resultsgraph = dcc.Graph(
                figure=fig,
                animate=False,
                # config={"scrollZoom": True},
                # style={"height": "60rem", "width": "100%"},
            )

            return [
                metadata_table,
                stats_table,
                paths_table,
                diff_table,
                resultsgraph,
                None,
            ]

    with server_is_started:
        server_is_started.notify()

    app.run_server(debug=False, port=port, host=host)
