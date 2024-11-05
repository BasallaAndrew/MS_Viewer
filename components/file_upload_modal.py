from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
from components import uploader

def render(app: Dash) -> html.Div:
    @app.callback(
        Output("upload-modal", "is_open"),
        Input('upload_button',"n_clicks"),
        State('upload-modal','is_open')
    )
    def toggle_upload_modal(n1,is_open):
        if n1:
            return not is_open
    return html.Div([
        dbc.Button("File Upload", id="upload_button", n_clicks=0),
        dbc.Modal([
            dbc.ModalHeader("File Upload"),
            dbc.ModalBody([
                html.Div(
                    children=[
                        uploader.render(app)
                    ]
                )
            ])
        ],
        id="upload-modal",
        size="xl",
        style={"max-width": "none"},
        is_open=False,
        className="xxl_modal",
        )
    ])