#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as ply
import numpy as np
import pandas as pd

def probability_graphic(file):
    # Gera um gráfico com a probabilidade de falha de propagação

    z_data = pd.read_csv(file)
    #np.genfromtxt(file,delimiter=',')

    P = z_data.as_matrix()
    #print(P)


    data = [
        go.Surface(
            z = P
        )
    ]
    layout = go.Layout(
        title='Probabilidade de falha na propagação do AP',
        autosize=False,
        width=600,
        height=600,
        margin=dict(
            l=10,
            r=10,
            b=10,
            t=10
        ),
        scene={"xaxis": {'title':"Vj [mV]", "tickfont": {"size": 10}, 'type': "linear"},
                    "yaxis": {'title': "N", "tickfont": {"size": 10},
                                "tickangle": 1},
                    "zaxis": {'title': "Pf",
                                "tickfont": {"size": 10}},
                    "camera": {"eye": {"x": 2, "y": 1, "z": 1.25}},
                    "aspectmode": "cube",
                    }
    )
    fig = go.Figure(data=data, layout=layout)
    ply.plot(fig, filename='Probabilidade_de_Falha_na_Propagacao.html')
    """
    trace = go.Scatter3d(
        z=P,
        marker=dict(
            size=4,
            color=P,
            colorscale='Viridis',
        ),
        line=dict(
            color='#1f77b4',
            width=1
        )
    )

    data = [trace]

    layout = dict(
        width=800,
        height=700,
        autosize=False,
        title='Failure Conductance Probability',
        scene=dict(
            xaxis=dict(
                title='Vj(mV)',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                title='N',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                title='Pf',
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            camera=dict(
                up=dict(
                    x=0,
                    y=0,
                    z=1
                ),
                eye=dict(
                    x=-1.7428,
                    y=1.0707,
                    z=0.7100,
                )
            ),
            aspectratio = dict( x=1, y=1, z=0.7 ),
            aspectmode = 'manual'
        ),
    )

    fig = dict(data=data, layout=layout)

    ply.plot(fig, filename='Probabilidade_de_Falha_na_Propagacao.html')
    """
