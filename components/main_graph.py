from dash import html, dcc, Output, Input
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd

# from .functions import generate_TIC

# from data.df import df

# retention_times, intensities, ms1_dict, ms2_dict = generate_TIC(r'data/022824_A_BSA_01.mzML')

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

        
        # figure = px.line(x= retention_times, y = intensities)
        # if clickData != None:
        #     figure.add_vline(x=clickData["points"][0]['x'],line_width=1, line_dash="dash", line_color="red")

        return 
    
    # @app.callback(
    #     [Output("align-graph", "figure")],  
    #     [
    #         Input("align-graph", "clickData")
    #     ],  
    # )
    # def update_main_graph(clickData):
    #     figure = px.line(x= retention_times, y = intensities)
    #     figure.add_vline(x=clickData["points"][0]['x'], line_width=1, line_dash="dash", line_color="red")
    #     return figure
    
    return html.Div([dcc.Graph(id="align-graph", responsive=True)], className="main")
