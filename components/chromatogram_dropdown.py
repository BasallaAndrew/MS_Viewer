from dash import Dash, html, dcc, Output, Input

def render (app: Dash) -> html.Div:
    chromatograms = ["TIC", "BPC", "XIC"]

    return html.Div(
        children=[
            html.H6("Chromatogram type"),
            dcc.Dropdown(
                id = "chromatogram-dropdown",
                options = [{"label": chrom, "value": chrom} for chrom in chromatograms],
                value='TIC',
                multi=False,
            ),
            html.Div(id='mz-range-container', children=[
                html.Label("Enter m/z range"),
                dcc.Input(id='mz-range', type='text', placeholder='e.g., 500-505',style={'width': '125px'})
            ])         
        ]
    )