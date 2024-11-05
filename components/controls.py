from dash import html, dcc
import dash_bootstrap_components as dbc
from components import file_upload_modal

def render(app) -> list:
    return html.Div(
        [
            html.H5("Interactive plotting"),
            html.P("Click on trace to extract data"),
            html.Hr(),
            html.Div(
                [
                    html.H5("Commands"),
                    dbc.Button(
                        "Display Graph",
                        id="display-graph",
                        color="primary",
                        style={"margin": "5px"},
                        n_clicks_timestamp="0",
                    ),
                ]
            ),
            file_upload_modal.render(app)
            
        ],
        className="controls",
    )
