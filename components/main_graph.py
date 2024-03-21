from dash import html, dcc, Output, Input
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd

from data.df import df


def render(app) -> html.Div:
    @app.callback(
        [Output("align-graph", "figure")],
        [
            Input("display-graph", "n_clicks"),
        ],
    )
    def display_main_graph(
        n_clicks,
    ):
        if not n_clicks:
            raise PreventUpdate

        figure = px.line(df, x=df.columns[0], y=df.columns[1:], template="plotly_white")

        figure.update_layout(
            autosize=True,
            width=1000,
            height=700,
            margin=dict(
                l=0,
                r=0,
                b=0,
                t=0,
            ),
        )
        return [figure]

    return html.Div([dcc.Graph(id="align-graph", responsive=True)], className="main")
