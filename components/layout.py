from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import pandas as pd

from components import main_graph, header, controls, outputs


def create_layout(app: Dash) -> html.Div:
    return html.Div(
        [
            header.render(app),
            controls.render(app),
            main_graph.render(app),
            outputs.render(app),
        ],
        className="grid-container",
    )
