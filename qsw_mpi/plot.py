#   QSW_MPI -  A package for parallel Quantum Stochastic Walk simulation.
#   Copyright (C) 2019 Edric Matwiejew
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import qsw_mpi.io as io
import qsw_mpi.measure as measure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.colors as colorz
from matplotlib.ticker import MaxNLocator

"""
Basic visualization of qsw_mpi time-series data and augmented graph structure.
"""

def population_lines(
        ax,
        pops,
        qsw_times,
        plot_times = [None, None],
        labels = False):
    """
    Line plot the vertex populations from a QSW time series.

    :param ax: Matplotlib :ref:`axis <projection-anchor>`

    :param pops: Site populations of a series as given by :meth:`~qsw_mpi.measure.populations`.
    :type pops: :math:`\\tilde{N}`, float, array

    :param qsw_times: Start, :math:`t_1`, and end, :math:`t_q`, times of the series.
    :type qsw_times: [float, float]

    :param plot_times: Start and end time for the plotting window.
    :type plot_times: [float, float], optional

    :param labels: If true, a legend labeling the vertices is included.
    :type labels: boolean, optional
    """
    steps = pops.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    plot_step_min = int((plot_times[0] - qsw_times[0])/h)
    plot_step_max = int((plot_times[1] - qsw_times[0])/h)

    ts = np.arange(qsw_times[0] + plot_step_min*h, qsw_times[0] + (plot_step_max)*h, h)

    for i in range(pops.shape[1]):
        if labels:
            ax.plot(ts, pops[plot_step_min:plot_step_max,i], label = str(i))
        else:
            ax.plot(ts, pops[plot_step_min:plot_step_max,i])

    if labels:
        ax.legend(title="vertex")

    plt.xlabel('time')
    plt.ylabel('population')

def coherence_lines(
        ax,
        node_pairs,
        cohs,
        qsw_times,
        plot_times = [None, None],
        labels = False):
    """
    Line plot the inter-vertex coherences from a QSW time series.

    :param node_pairs: Array labeling the node-pairs as given by :meth:`~qsw_mpi.measure.coherences`.
    :type node_pairs: integer

    :param cohs: Inter-vertex coherences as given by :meth:`~qsw_mpi.measure.coherences`.
    :type cohs: integer

    :param qsw_times: Start, :math:`t_1`, and end, :math:`t_q`, times of the series.
    :type qsw_times: [float, float]

    :param plot_times: Start and end time for the plotting window.
    :type plot_times: [float, float], optional

    :param labels: If true, a legend labeling the vertex pairs is included.
    :type labels: boolean, optional
    """
    steps = cohs.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    plot_step_min = int((plot_times[0] - qsw_times[0])/h)
    plot_step_max = int((plot_times[1] - qsw_times[0])/h)

    ts = np.arange(qsw_times[0] + plot_step_min*h, qsw_times[0] + plot_step_max*h, h)

    for i in range(cohs.shape[1]):
        if labels:
            ax.plot(ts, cohs[plot_step_min:plot_step_max,i], label = str((node_pairs[0][i], node_pairs[1][i])))
        else:
            ax.plot(ts, cohs[plot_step_min:plot_step_max,i])

    if labels:
        ax.legend(title = "vertex pairs")

    plt.ylabel(r'coherence')
    plt.xlabel('time')

def population_bars(
        ax,
        pops,
        qsw_times,
        plot_times = [None, None]):

    """
    3D bar plot the vertex populations from a QSW time series.

    :param ax: Matplotlib axis with a :ref:`3D projection <projection-anchor>`

    :param pops: Site populations of a series as given by :meth:`~qsw_mpi.measure.populations`.
    :type pops: :math:`\\tilde{N}`, float, array

    :param qsw_times: Start, :math:`t_0`, and end, :math:`t_q`, times of the series.
    :type qsw_times: [float, float]

    :param plot_times: Start and end time for the plotting window.
    :type plot_times: [float, float], optional
    """
    steps = pops.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    x = np.arange(0,pops.shape[1],1)

    h = (qsw_times[1] - qsw_times[0])/float(steps)
    y = np.full(pops.shape[1],0)

    it = iter(range(steps))
    next(it,None)

    for i in it:
        x = np.vstack((x,np.arange(0,pops.shape[1],1)))
        y = np.vstack((y, np.full(pops.shape[1], i)))

    x = np.ndarray.flatten(x)
    y = h*np.ndarray.flatten(y)

    z = np.zeros(x.shape)

    dx = np.full(x.shape, 0.2)
    dy = np.full(x.shape, 0.2)
    dz = np.ndarray.flatten(pops)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_xlabel('vertex')
    ax.set_ylabel('time')
    ax.set_zlabel('population')

    fracs = dz.astype(float)/dz.max()
    norm = colorz.Normalize(fracs.min(), fracs.max())
    colors = cm.jet(norm(fracs))

    ax.bar3d(x, y, z, dx, dy, dz, color=colors)

    ax.view_init(30, 115)

def coherence_bars(
        ax,
        node_pairs,
        cohs,
        qsw_times,
        plot_times = [None, None]):
    """
    3D bar plot the vertex populations from a QSW time series.

    :param ax: Matplotlib axis with a :ref:`3D projection <projection-anchor>`

    :param node_pairs: Array labeling the node-pairs as given by :meth:`~qsw_mpi.measure.coherences`.
    :type node_pairs: integer

    :param cohs: Inter-vertex coherences as given by :meth:`~qsw_mpi.measure.coherences`.
    :type cohs: integer

    :param qsw_times: Start, :math:`t_1`, and end, :math:`t_q`, times of the series.
    :type qsw_times: [float, float]

    :param plot_times: Start and end time for the plotting window.
    :type plot_times: [float, float], optional

    """

    steps = cohs.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    x = np.arange(0,node_pairs[0].shape[0],1)
    y = np.full(cohs.shape[1],0)

    it = iter(range(steps))
    next(it,None)

    for i in it:
        x = np.vstack((x,np.arange(0,node_pairs[0].shape[0],1)))
        y = np.vstack((y, np.full(cohs.shape[1], i)))

    x = np.ndarray.flatten(x)
    y = h*np.ndarray.flatten(y)

    z = np.zeros(x.shape)

    dx = np.full(x.shape, 0.2)
    dy = np.full(x.shape, 0.2)
    dz = np.ndarray.flatten(cohs)

    plt.xticks(x, [str(node_pair) for node_pair in zip(node_pairs[0], node_pairs[1])], rotation = 45, ha = 'right')

    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(-7.7)

    ax.set_xlabel('vertex pairs')
    ax.set_ylabel('time')
    ax.set_zlabel('coherence')

    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colorz.Normalize(fracs.min(), fracs.max())
    colors = cm.jet(norm(fracs))

    ax.bar3d(x, y, z, dx, dy, dz, color=colors)

    ax.view_init(30, 115)

def graph(
        G,
        sources = None,
        sinks = None,
        layout = None,
        graph_color = 'yellow',
        source_color = 'lightgreen',
        sink_color = 'pink',
        node_labels = True,
        node_font_size = 8,
        node_size = 400,
        legend = True,
        legend_font_size = 12,
        legend_key_size = 300,
        legend_label_graph = 'graph',
        legend_label_source = 'source',
        legend_label_sink = 'sink'):

    """
    Visualize a graph with attached sources and sinks.

    :param G: Networkx graph class.

    :param sources: Source array tuple as defined in :meth:`~qsw_mpi.MPI.walk`.
    :type sources: ([M], [M]), (integer, float), tuple, optional

    :param sinks: Sink array tuple as defined in :meth:`~qsw_mpi.MPI.walk`.
    :type sinks: ([M], [M]), (integer, float), tuple, optional

    :param layout: Networkx graph layout.
    :type layout: optional

    :param graph_color: Colour of the graph vertices.
    :type graph_color: string, optional

    :param source_color: Colour of the source vertices.
    :type source_color: string, optional

    :param sink_color: Colour of the sink vertices.
    :type sink_color: string, optional

    :param node_labels: If True label the vertices.
    :type node_labels: boolean, optional

    :param node_font_size: Font size of the node labels.
    :type node_font_size: integer, optional

    :param node_size: Vertex size.
    :type node_size: integer, optional

    :param legend: If True include legend defining the graph, source and sink vertices.
    :type legend: boolean, optional

    :param legend_font_size: Font size of the legend text.
    :type legend_font_size: integer, optional

    :param legend_key_size: Size of the legend markers.
    :type legend_key_size: integer, optional

    :param legend_label_graph: Name for the graph vertices.
    :type legend_label_graph: string, optional

    :param legend_label_source: Name for the source vertices.
    :type legend_label_source: string, optional

    :param legend_label_sink: Name fot the sink vertices.
    :type legend_label_sink: string, optional
    """

    def create_legend():

        plt.rc('legend', fontsize = legend_font_size)
        lgnd = plt.legend(scatterpoints=1)
        lgnd.legendHandles[0]._sizes = [legend_key_size]

        if sources is not None:
            lgnd.legendHandles[1]._sizes = [legend_key_size]
        if sinks is not None:
            lgnd.legendHandles[2]._sizes = [legend_key_size]

    if (sinks is None) and (sources is None):
        Gaug = G
    else:
        Gaug = G.copy()

    nx.set_node_attributes(Gaug, values=" ", name='genus')

    if sources is not None:
        n = nx.number_of_nodes(Gaug)
        for i, source in enumerate(sources[0]):
            sourceIdx = n + i
            Gaug.add_node(sourceIdx,**{'genus':'source'})
            Gaug.add_edge(sourceIdx, source)

    if sinks is not None:
        n = nx.number_of_nodes(Gaug)
        for i, sink in enumerate(sinks[0]):
            sinkIdx = n + i
            Gaug.add_node(sinkIdx, **{'genus':'sink'})
            Gaug.add_edge(sinkIdx, sink)

    graph_node_list = []
    graph_color_map = []
    graph_labels = {}

    sink_node_list = []
    sink_color_map = []
    sink_labels = {}

    source_node_list = []
    source_color_map = []
    source_labels = {}

    for node in Gaug:
        if Gaug.nodes[node]['genus'] is 'sink':
            sink_labels[node] = node
            sink_color_map.append(sink_color)
            sink_node_list.append(node)
        elif Gaug.nodes[node]['genus'] is 'source':
            source_labels[node] = node
            source_color_map.append(source_color)
            source_node_list.append(node)
        else:
            graph_labels[node] = node
            graph_color_map.append(graph_color)
            graph_node_list.append(node)

    if layout is None:
        layout = nx.spring_layout(Gaug, iterations = 50)

    plt.axis('off')

    nx.draw_networkx_nodes(
            Gaug,
            pos = layout,
            node_color = graph_color_map,
            nodelist = graph_node_list,
            label = legend_label_graph,
            node_shape = 'o',
            node_size = node_size)

    if node_labels:
        nx.draw_networkx_labels(
                Gaug,
                pos = layout,
                nodelist = graph_node_list,
                font_size = node_font_size)

    nx.draw_networkx_edges(
            Gaug,
            pos = layout,
            nodelist = graph_node_list)

    if sources is not None:
        nx.draw_networkx_nodes(
                Gaug,
                pos = layout,
                node_color = source_color_map,
                nodelist = source_node_list,
                label = legend_label_source,
                node_shape = '*',
                node_size = node_size)

        nx.draw_networkx_edges(
                Gaug,
                pos = layout,
                nodelist = source_node_list)

    if sinks is not None:
        nx.draw_networkx_nodes(
                Gaug,
                pos = layout,
                node_color = sink_color_map,
                nodelist = sink_node_list,
                label = legend_label_sink,
                node_shape = 'H',
                node_size = node_size)

        nx.draw_networkx_edges(
                Gaug,
                pos = layout,
                nodelist = sink_node_list)

    if legend:
        create_legend()

    elif (sources is not None) or (sinks is not None):
        create_legend()

def animate(
        fig,
        G,
        populations,
        start,
        end,
        filename = 'animation',
        save = True,
        node_size = 3000,
        framerate = 25,
        graph_color = 'orange',
        layout = None,
        title = None,
        sources = None,
        sinks = None,
        use_labels = False):

    """
    Create an animation of a QSW.

    :param fig: Matplotlib figure object.

    :param G: Networkx graph object.

    :param populations: Vertex populations as given by :meth:`~qsw_mpi.measure.populations`.
    :type populations: :math:`\\tilde{N}`, float, array

    :param start: Start time of the walk.
    :type start: float

    :param end: End time of the walk.
    :type end: float

    :param filename: Output filename.
    :type filename: string

    :param save: If True save the animation to disk.
    :type save: boolean, optional

    :param node_size: Maximum vertex size.
    :type node_size: integer, optional

    :param framerate: Animation framerate.
    :type framerate: integer, optional

    :param graph_color: Colour of the graph vertices.
    :type graph_color: string, optional

    :param layout: Networkx graph layout.
    :type layout: optional

    :param title: Animation title.
    :type title: string, optional

    :param sources: Source array tuple as defined in :meth:`~qsw_mpi.MPI.walk`.
    :type sources: ([M], [M]), (integer, float), tuple, optional

    :param sinks: Sink array tuple as defined in :meth:`~qsw_mpi.MPI.walk`.
    :type sinks: ([M], [M]), (integer, float), tuple, optional

    :param use_labels: If True label the vertices.
    :type use_labels: boolean, optional

    """

    if (sinks is None) and (sources is None):
        Gaug = G
    else:
        Gaug = G.copy()

    labels = {}
    color_map = []

    nx.set_node_attributes(Gaug, name='genus', values="")

    if sources is not None:
        n = nx.number_of_nodes(G)
        for i, source in enumerate(sources[0]):
            sourceIdx = n + i
            Gaug.add_node(sourceIdx,**{'genus':'source'})
            Gaug.add_edge(sourceIdx, source)

    if sinks is not None:
        n = nx.number_of_nodes(Gaug)
        for i, sink in enumerate(sinks[0]):
            sinkIdx = n + i
            Gaug.add_node(sinkIdx, **{'genus':'sink'})
            Gaug.add_edge(sinkIdx, sink)

    node_list = []
    for node in Gaug:
        if Gaug.nodes[node]['genus'] is 'sink':
            labels[node] = 'sink ' + str(node)
            color_map.append('red')
        elif Gaug.nodes[node]['genus'] is 'source':
            labels[node] = 'source ' + str(node)
            color_map.append('green')
        else:
            labels[node] = node
            color_map.append(graph_color)
        node_list.append(node)

    if layout is None:
        layout = nx.spring_layout(Gaug, iterations = 1000)

    if title is None:
        title = ""

    steps = populations.shape[0]
    deltat = float(end)/steps

    nodesize = np.full(Gaug.number_of_nodes(), 0)

    plt.axis('off')
    plt.title(title)
    plt.tight_layout()

    nx.draw_networkx_edges(
            Gaug,
            pos=layout,
            node_list = node_list,
            node_size = nodesize)

    if use_labels:
        nx.draw_networkx_labels(
                Gaug,
                pos=layout,
                node_list = node_list,
                labels = labels,
                node_size = nodesize)

    ann_list = []
    def animate(i):

        for a in ann_list:
            a.remove()
        ann_list[:] = []

        t= start + np.mod(i,steps)*deltat
        time = plt.annotate("t = " + str(round(t,2)), [20,20],\
        xycoords='figure pixels')
        ann_list.append(time)

        for j, node in enumerate(Gaug):
            nodesize[j] = node_size*populations[np.mod(i,steps)][node]

        nx.draw_networkx_nodes(
                Gaug,
                pos=layout,
                node_size = nodesize,
                node_color=color_map,
                node_list = node_list)

    miliseconds= math.ceil(1000/framerate)

    ani = animation.FuncAnimation(
            fig,
            animate,
            interval=miliseconds,
            frames = steps)

    if save:
        ani.save(filename + '.gif', fps=framerate, writer='imagemagick')
        plt.clf()
    else:
        plt.show()

