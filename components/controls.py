from dash import html, dcc
import dash_bootstrap_components as dbc


def render() -> list:
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
        ],
        className="controls",
    )
