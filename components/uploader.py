from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_bootstrap_components as dbc

def render(app: Dash) -> html.Div:
    return html.Div(
        children= [dcc.Upload(id='upload-data',
                              children=html.Div([
                                  'Drag and Drop or ',
                                  html.A('Select Files')
                              ]),
                              style={
                                  'width':'100%',
                                  'height':'60px',
                                  'lineHeight': '60px',
                                  'borderWidth' : '1px',
                                  'borderStyle': 'dashed',
                                  'borderRadius': '5px',
                                  'textAlign' : 'center',
                                  'margin': '10px'
                              },
                              multiple=False)

        ]
    )