import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import freeqsw.io as io
import freeqsw.measure as measure

def graph(  G,\
            filename, \
            sources = None,\
            sinks = None, \
            title = None, \
            title_font_size = 16, \
            layout = None,\
            graph_color = 'yellow', \
            source_color = 'lightgreen', \
            sink_color = 'pink', \
            size = (5,5), \
            format = 'png', \
            node_labels = True, \
            node_font_size = 8, \
            node_size = 400, \
            legend = None, \
            legend_font_size = 12, \
            legend_key_size = 300, \
            legend_label_graph = 'graph', \
            legend_label_source = 'source', \
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
        #fixed_positions = {0:(0,0)}
        #fixed_nodes = fixed_positions.keys()
        #layout = nx.spring_layout(Gaug, iterations = 50, fixed = fixed_nodes, pos=fixed_positions)
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

    plt.savefig(filename + '.' + format, format = format)
    plt.clf()


def animate(G,
            populations,
            start,
            end,
            filename, \
            node_size = 1000,\
            framerate = 25,\
            graph_color = 'orange',\
            layout = None,\
            title = None,\
            sources = None,\
            sinks = None, \
            size = (5,5), \
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
        #fixed_positions = {0:(0,0)}
        #fixed_nodes = fixed_positions.keys()
        #layout = nx.spring_layout(Gaug, iterations = 500, fixed = fixed_nodes, pos=fixed_positions)
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
