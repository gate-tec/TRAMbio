import plotly.graph_objects as go
from nicegui import ui

def init_plot():
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 10], y=[1, 1],
                             line=dict(color="rgb(189,189,189)", width=2),
                             showlegend=False,
                             mode='lines', hoverinfo='none'))
    fig.add_trace(go.Scatter(x=[0, 10, 10, 0], y=[1.3, 1.3, 0.7, 0.7],
                             fill='toself', fillcolor="rgba(255,255,255,0)",
                             line_color="rgba(255,255,255,0)", mode='lines',
                             showlegend=False, hoverinfo='none'))

    fig.add_trace(go.Scatter(x=[5,5,8,8], y=[0.7,1.3,1.3,0.7],
                        fill='toself', fillcolor='darkviolet',
                        hoveron = 'fills', # select where hover is active
                        line_color='darkviolet', mode='lines',
                        text="Points + Fills",
                        hoverinfo='text+x+y'))
    fig.add_shape(
        type="rect",
        fillcolor='LightGreen',
        x0=0,
        y0=0,
        x1=2,
        y1=2,
        label=dict(text="Text at 45", textangle=45),
    )

    fig.add_shape(
        type="rect",
        fillcolor='Gold',
        x0=3,
        y0=0,
        x1=5,
        y1=2,
        label=dict(text="Text at -45", textangle=-45),
    )

    plot = ui.plotly(fig).classes("w-full")
    plot.on("plotly_hover", lambda x: print("Hover:", x))