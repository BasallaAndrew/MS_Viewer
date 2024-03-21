import dash_bootstrap_components as dbc
from dash import Dash

from components.layout import create_layout

app = Dash(external_stylesheets=[dbc.themes.MINTY])
app.title = "Simple Graph"
server = app.server

app.layout = create_layout(app)

if __name__ == "__main__":
    app.run_server(debug=True)
