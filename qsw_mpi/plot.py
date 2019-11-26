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

def population_lines(pops, qsw_times, plot_times = [None, None], labels = False, figsize = (5,4)):

    fig = plt.figure(figsize=figsize)

    steps = pops.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    plot_step_min = int((plot_times[0] - qsw_times[0])/h)
    plot_step_max = int((plot_times[1] - qsw_times[0])/h - 1)

    ts = np.arange(qsw_times[0] + plot_step_min*h, qsw_times[0] + (plot_step_max)*h, h)

    for i in range(pops.shape[1]):
        if labels:
            plt.plot(ts, pops[plot_step_min:plot_step_max,i], label = str(i))
        else:
            plt.plot(ts, pops[plot_step_min:plot_step_max,i])

    if labels:
        plt.legend(title="vertex")
    plt.xticks([0,4,8,12,16,20])
    plt.ylabel(r'population',fontsize = 14)
    plt.xlabel('t',fontsize = 14)

    return fig

def coherence_lines(node_pairs, cohs, qsw_times, plot_times = [None, None], labels = False, figsize = (5,4)):

    fig = plt.figure(figsize=figsize)

    steps = cohs.shape[0]

    h = (qsw_times[1] - qsw_times[0])/float(steps)

    if plot_times[0] is None:
        plot_times[0] = qsw_times[0]
    if plot_times[1] is None:
        plot_times[1] = qsw_times[1]

    plot_step_min = int((plot_times[0] - qsw_times[0])/h)
    plot_step_max = int((plot_times[1] - qsw_times[0])/h - 1)

    ts = np.arange(qsw_times[0] + plot_step_min*h, qsw_times[0] + plot_step_max*h, h)

    for i in range(cohs.shape[1]):
        if labels:
            plt.plot(ts, cohs[plot_step_min:plot_step_max,i], label = str((node_pairs[0][i], node_pairs[1][i])))
        else:
            plt.plot(ts, cohs[plot_step_min:plot_step_max,i])

    if labels:
        plt.legend(title = "vertex pairs")
    plt.xticks([0,4,8,12,16,20])
    plt.ylabel(r'coherence',fontsize = 14)
    plt.xlabel('t',fontsize = 14)

    return fig

def population_bars(pops, t1, t2, t_tick_freq = None, t_round = 2, figsize = (5,4)):

    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d')

    steps = pops.shape[0]

    pops = np.flip(pops)

    x = np.arange(0,pops.shape[1],1)

    h = (t2 - t1)/float(steps)
    y = np.full(pops.shape[1],0)

    it = iter(range(steps))
    next(it,None)

    for i in it:
        x = np.vstack((x,np.arange(0,pops.shape[1],1)))
        y = np.vstack((y, np.full(pops.shape[1], i)))

    x = np.ndarray.flatten(x)
    y = np.ndarray.flatten(y)


    z = np.zeros(x.shape)

    dx = np.full(x.shape, 0.2)
    dy = np.full(x.shape, 0.2)
    dz = np.ndarray.flatten(pops)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if t_tick_freq is not None:
        plt.yticks([i for i in range(0,steps,t_tick_freq)], [str(round(t,t_round)) for t in np.arange(h,steps,t_tick_freq*h)])
    else:
        plt.yticks([i for i in range(0,steps)], [str(round(t,2)) for t in np.arange(h,steps,h)])
    plt.yticks([0,40,80,120,160,200],[20,16,12,8,4,0])
    plt.ylim(0,steps + 1)
    plt.xticks([0,1,2,3,4,5],[5,4,3,2,1,0])
    ax.set_xlabel('vertex',fontsize = 14)
    ax.set_ylabel('t')
    ax.set_zlabel(r'population',fontsize = 14)
    fracs = dz.astype(float)/dz.max()
    norm = colorz.Normalize(fracs.min(), fracs.max())
    colors = cm.jet(norm(fracs))

    ax.bar3d(x, y, z, dx, dy, dz, color=colors)

    return fig

def coherence_bars(node_pairs, cohs, t1, t2, t_tick_freq = None, t_round = 2, figsize = (5,4)):

    fig = plt.figure(figsize=figsize)
    ax = fig.gca(projection='3d')

    steps = cohs.shape[0]

    x = np.arange(0,node_pairs[0].shape[0],1)
    h = (t2 - t1)/float(steps)
    y = np.full(cohs.shape[1],0)

    it = iter(range(steps))
    next(it,None)

    for i in it:
        x = np.vstack((x,np.arange(0,node_pairs[0].shape[0],1)))
        y = np.vstack((y, np.full(cohs.shape[1], i)))

    x = np.ndarray.flatten(x)
    y = np.ndarray.flatten(y)

    z = np.zeros(x.shape)

    dx = np.full(x.shape, 0.2)
    dy = np.full(x.shape, 0.2)
    dz = np.ndarray.flatten(cohs)

    plt.xticks(x, [str(node_pair) for node_pair in zip(node_pairs[0], node_pairs[1])], rotation = 45, ha = 'right')

    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(-7.7)

    if t_tick_freq is not None:
        plt.yticks([i for i in range(0,steps,t_tick_freq)], [str(round(t,t_round)) for t in np.arange(h,steps,t_tick_freq*h)])
    else:
        plt.yticks([i for i in range(0,steps)], [str(round(t,2)) for t in np.arange(h,steps,h)])
    plt.yticks([0,40,80,120,160,200],[0,4,8,12,16,20])
    plt.ylim(steps + 1, 0)

    ax.set_xlabel('vertex pairs',fontsize = 13)
    ax.set_ylabel('t')
    ax.set_zlabel(r'coherence',fontsize = 14)

    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colorz.Normalize(fracs.min(), fracs.max())
    colors = cm.jet(norm(fracs))

    ax.bar3d(x, y, z, dx, dy, dz, color=colors)

    return fig


def graph(
        G,
        sources = None,
        sinks = None,
        title = None,
        title_font_size = 16,
        layout = None,
        graph_color = 'yellow',
        source_color = 'lightgreen',
        sink_color = 'pink',
        size = (5,5),
        format = 'png',
        node_labels = True,
        node_font_size = 8,
        node_size = 400,
        legend = None,
        legend_font_size = 12,
        legend_key_size = 300,
        legend_label_graph = 'graph',
        legend_label_source = 'source',
        legend_label_sink = 'sink'):

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

    nx.set_node_attributes(Gaug, name='genus', values="")

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
        if Gaug.node[node]['genus'] is 'sink':
            sink_labels[node] = node
            sink_color_map.append(sink_color)
            sink_node_list.append(node)
        elif Gaug.node[node]['genus'] is 'source':
            source_labels[node] = node
            source_color_map.append(source_color)
            source_node_list.append(node)
        else:
            graph_labels[node] = node
            graph_color_map.append(graph_color)
            graph_node_list.append(node)

    if layout is None:
        layout = nx.spring_layout(Gaug, iterations = 50)

    fig = plt.figure(figsize = (size[0], size[1]), dpi = 80)

    plt.axis('off')

    if title is not None:
        plt.rc('figure', titlesize = title_font_size)
        plt.title(title)

    nx.draw_networkx_nodes(Gaug, pos = layout, node_color = graph_color_map, \
        nodelist = graph_node_list, label = legend_label_graph, node_shape = 'o', node_size = node_size)

    if node_labels:
        nx.draw_networkx_labels(Gaug, pos = layout, nodelist = graph_node_list, font_size = node_font_size)

    nx.draw_networkx_edges(Gaug, pos = layout, nodelist = graph_node_list)

    if sources is not None:
        nx.draw_networkx_nodes(Gaug, pos = layout, node_color = source_color_map, \
            nodelist = source_node_list, label = legend_label_source, node_shape = '*', node_size = node_size)

        nx.draw_networkx_edges(Gaug, pos = layout, nodelist = source_node_list )

    if sinks is not None:
        nx.draw_networkx_nodes(Gaug, pos = layout, node_color = sink_color_map, \
            nodelist = sink_node_list, label = legend_label_sink, node_shape = 'H', node_size = node_size)

        nx.draw_networkx_edges(Gaug, pos = layout, nodelist = sink_node_list)

    if legend is not None:

        if legend:
            create_legend()

    elif (sources is not None) or (sinks is not None):
        create_legend()

    return fig

def animate(
        G,
        populations,
        start,
        end,
        filename,
        node_size = 1000,
        framerate = 25,
        graph_color = 'orange',
        layout = None,
        title = None,
        sources = None,
        sinks = None,
        size = (5,5),
        use_labels = False):

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
        if Gaug.node[node]['genus'] is 'sink':
            labels[node] = 'sink ' + str(node)
            color_map.append('red')
        elif Gaug.node[node]['genus'] is 'source':
            labels[node] = 'source ' + str(node)
            color_map.append('green')
        else:
            labels[node] = node
            color_map.append(graph_color)
        node_list.append(node)

    if layout is None:
        layout = nx.spring_layout(Gaug, iterations = 500)

    if title is None:
        title = ""

    steps = populations.shape[0]
    deltat = float(end)/steps

    fig = plt.figure(figsize = size)

    nodesize = np.full(Gaug.number_of_nodes(), 0)

    plt.axis('off')
    plt.title(title)
    plt.tight_layout()

    nx.draw_networkx_edges(Gaug, pos=layout, node_list = node_list, node_size = nodesize)

    if use_labels:
        nx.draw_networkx_labels(Gaug, pos=layout, node_list = node_list, labels = labels, node_size = nodesize)

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

        nx.draw_networkx_nodes(Gaug,pos=layout,\
                node_size = nodesize, node_color=color_map, node_list = node_list)

    miliseconds= int(1000/framerate)
    ani = animation.FuncAnimation(fig, animate, interval=miliseconds, frames = steps)
    ani.save(str(filename) + '.gif', writer='imagemagick', fps=framerate, bitrate = -1)
    plt.clf()
