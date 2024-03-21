from dash import html, dcc, Input, Output
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import json

from data.df import df

styles = {"pre": {"border": "thin lightgrey solid", "overflowX": "scroll"}}


def render(app) -> html.Div:
    @app.callback(
        Output("click-data", "children"),
        Output("trace-number-clicked", "children"),
        Input("align-graph", "clickData"),
    )
    def display_click_data(clickData):
        if clickData and clickData["points"]:
            trace_index = int(clickData["points"][0]["curveNumber"]) + 1
            trace_name_clicked = df.columns[trace_index]
            return json.dumps(clickData, indent=2), trace_name_clicked
        raise PreventUpdate

    return html.Div(
        [
            html.H5("Outputs"),
            dcc.Markdown(
                """
                **Click Data**

                Click on points in the graph.
            """
            ),
            html.Pre(id="click-data", style=styles["pre"]),
            html.Pre(id="trace-number-clicked"),
        ],
        className="outputs",
    )
