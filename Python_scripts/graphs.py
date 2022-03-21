import cobra
import copy
import matplotlib.pyplot as plt

from upsetplot import from_memberships, plot

from utils import list_reactions_cobra, cobra_compatibility, find_files, write_file


def get_clusters(cluster_list):
    """Function to create every individual cluster depending on
    the number of values given to the UpSetPlot function.

    PARAMS:
        cluster_list (list) -- the list of subjects for the upSetPlot.

    RETURNS:
        final_res (list of lists) -- a list of all the possible clusters.
    """

    res = []
    final_res = []
    for i in range(len(cluster_list) - 1):
        if i == 0:
            for x in cluster_list:
                z = cluster_list.index(x)
                for j in range(len(cluster_list) - z - 1):
                    res.append([x, cluster_list[z + j + 1]])
            [final_res.append(k) for k in res]
        else:
            res = get_sub_clusters(cluster_list, res)
            [final_res.append(r) for r in res]
    return final_res


def get_sub_clusters(cluster_list, res):
    """Sub-function of the clusters (algorithmic architecture)."""

    sub_res = []
    for y in res:
        z = cluster_list.index(y[len(y) - 1])
        for i in range(z + 1, len(cluster_list)):
            x = copy.deepcopy(y)
            x.append(cluster_list[i])
            sub_res.append(x)
    return sub_res


def similarity_count(data, args, others):
    """Function which is part of the process to make the UpSetPlot,
    counting and returning the similarities between the clusters."""

    cluster_set = set.intersection(*args)
    for i in others:
        cluster_set = cluster_set.difference(set(data[i]))
    return cluster_set, len(cluster_set)


def make_upsetplot(directory, name, data, title, remove_zero=True, show_plot=False):
    """Function to make an UpSetPlot.
    Need this three other functions : similarity_count(), get_clusters(), get_sub_clusters().

    PARAMS:
        directory (str) -- the directory to save the result.
        name (str) -- name of the file to save.
        data -- the dictionary containing the organisms as keys
        and the genes/reactions/others to treat for the UpSetPlot.
        title (str) -- title of the graph.
        remove_zero (boolean) -- False if you want to plot the 0 values intersections.
        show_plot (boolean) -- True if you want a graph to pop.
    """

    clusters = get_clusters(list(data.keys()))
    [clusters.insert(0, [key]) for key in data.keys()]
    count = []
    log = ""
    for c in clusters:
        others = list(data.keys())
        list_inter = []
        for x in c:
            others.remove(x)
            list_inter.append(set(data[x]))
        cluster_data, sim_count = similarity_count(data, list_inter, others)
        count.append(sim_count)
        for i in c:
            log += i + " "
        log += " (" + str(sim_count) + ") :\n\n"
        for i in cluster_data:
            log += cobra_compatibility(str(i)) + "\n"
        log += "\n------\n\n"
    write_file(directory + name + ".log", log, False)
    # removing 0 values
    if remove_zero:
        while 0 in count:
            x = count.index(0)
            count.pop(x)
            clusters.pop(x)
    my_upsetplot = from_memberships(clusters, count)
    plot(my_upsetplot, show_counts='%d', totals_plot_elements=3)
    plt.suptitle(title)
    plt.savefig(directory + name + ".pdf")
    if show_plot:
        plt.show()


if __name__ == '__main__':
    # Example of use, put all the .json networks you want to compare in a directory and point to this directory :
    wd = '/home/asa/Bureau/Graphs/'
    dicoUpset = {}
    for network in find_files(wd, '.json'):
        model = cobra.io.load_json_model(wd + network)
        dicoUpset[model.id] = list_reactions_cobra(model)

    # Choose a name, an explanation and let's go !
    make_upsetplot(wd, "Name", dicoUpset,  "Explanation")
