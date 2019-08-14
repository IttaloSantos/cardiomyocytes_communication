#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as ply

# Gera um gráfico com a probabilidade de falha de propagação

trace = np.arange(N)
for i in range(N):
    trace[i] = go.Scatter3d(
        x=x[:, i], y=y, z=z,
        marker=dict(
            size=4,
            color=z,
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
