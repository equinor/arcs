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

def keys_by_depth(dict_, depth=0, output=None):
    if output is None:
        output = {}
    if depth not in output:
        output[depth] = set()
    for key in dict_:
        output[depth].add(key)
        if isinstance(dict_[key], dict):
            keys_by_depth(dict_[key], depth + 1, output)
    return(output)
    

def _markdown_compound(_string):
    md = []
    for i in _string:
        try:
            int(i)
            md.append("<sub>{}</sub>".format(int(i)))
        except Exception:
            md.append(i)
    return ("".join(md))
    

def format_concentrations_table_for_dash(dataframe):
    dataframe = dataframe.T
    dataframe = dataframe.T.rename(
        {x: _markdown_compound(x) for x in dataframe}
    )
    dataframe = dataframe.reset_index()
    dataframe = dataframe.rename(
        {
            "index": "compound",
            "initial": "initial (ppm)",
            "mean": "mean (ppm)",
            "diff": "change (ppm)",
            "sem": "sem (ppm)",
            "std": "std (ppm)",
            "var": "var (ppm)",
        }
    )
    formatted_table = dash_table.DataTable(
        columns=[
            {"name": i,
                "id": i,
                "type": "text",
                "presentation": "markdown"}
            for i in dataframe.columns
        ],
        data=dataframe.to_dict("records"),
        style_as_list_view=False,
        cell_selectable=False,
        style_cell={
            "font_family": "helvetica",
            "align": "left",
            'padding-right': '10px',
            'padding-left': '10px',
            'text-align': 'left',
            'marginLeft': 'auto',
            'marginRight': 'auto'
        },
        style_header={
            "font_family": "helvetica",
            "align": "left",
            'padding-right': '10px',
            'padding-left': '10px',
            'text-align': 'left',
            'marginLeft': 'auto',
            'marginRight': 'auto'
        },
        style_table={
            "overflow": "scroll",
        },
        markdown_options={"html": True, "link_target": "_self"},
    )
    return(formatted_table)


def format_statistics_table_for_dash(dataframe):
    formatted_stats_table = dash_table.DataTable(
        columns=[
            {
                "name": "Reactions",
                        "id": "index",
                        "type": "text",
                        "presentation": "markdown",
            },
            #{
            #    "name": "k",
            #            "id": "k",
            #            "type": "text",
            #            "presentation": "markdown",
            #},
            {
                "name": "Frequency",
                        "id": "frequency",
                        "type": "text",
                        "presentation": "markdown",
            },
        ],
        data=dataframe.to_dict("records"),
        style_as_list_view=False,
        cell_selectable=False,
        style_cell={
            "font_family": "helvetica",
            "align": "left",
            'padding-right': '10px',
            'padding-left': '10px',
            'text-align': 'left',
            'marginLeft': 'auto',
            'marginRight': 'auto'
        },
        style_header={
            "font_family": "helvetica",
            "align": "left",
            'padding-right': '10px',
            'padding-left': '10px',
            'text-align': 'left',
            'marginLeft': 'auto',
            'marginRight': 'auto'
        },
        style_table={
            "overflow": "scroll",
        },
        # fixed_rows={"headers": True},
        markdown_options={"html": True, "link_target": "_self"},
    )
    return(formatted_stats_table)

def format_bar_chart(dataframe):

    #dataframe.reset_index(inplace=True)
    
    dataframe = pd.DataFrame(
        {
            "comps": list(dataframe.T.keys()),
            "diff": dataframe['diff'].values,
            "sem": dataframe['sem'].values,
        }
    )
    
    fig = px.bar(
        dataframe,
        x="comps",
        y="diff",
        error_y="sem",
        #error_y_minus="sem_minus",
        labels={"comps": "", "diff": "\u0394 ppm"},
        color="diff",
        color_continuous_scale="tropic_r",
        hover_data={
            "diff": False,
            "comps": False,
            #"variance": False,
            "error": (":.2f", dataframe["sem"]),
            "specie": dataframe["comps"],
            "PPM": (":.2f", dataframe["diff"]),
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
    fig.update_xaxes(
        showgrid=False, 
        tickangle=-60,
        tickmode="linear")
            #try:
            #    dtick=int(int(ymax - ymin) / 10)
            #except:
            #    dtick = None
            #    pass
            #fig.update_yaxes(
            #    showgrid=True,
            #    tickmode="linear",
            #    range=[ymin - 2, ymax + 2],
            #    dtick=dtick,
            #)
    
    resultsgraph=dcc.Graph(
                figure=fig,
                animate=False,
                #config={"scrollZoom": True},
                #style={"height": "60rem", "width": "100%"},
            )
    return(resultsgraph)
