from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
from .functions import parse_mzml_to_dict

def render(app: Dash) -> html.Div:
    
    @app.callback(
    [Output('stored-data', 'data'),  # Save parsed data to dcc.Store
        Output('file-status', 'children')],  # Display feedback
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
    )
    def parse_and_store_data(contents, filename):
        if contents is None:
            return None, "No file uploaded."

        data = parse_mzml_to_dict(filename)
        return data, f"File '{filename}' uploaded successfully."

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
                              multiple=False),
                              dcc.Store(id='stored-data'),
                              html.Div(id='file-status')

        ]
    )

