#!/usr/bin/env python2.7
"""
Cluster some graphs (using vg compare for pairwise distances) based on their similarity.
Current implementation : neighbour joining tree using Jaccard distance matrix
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, math
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import pylab
import networkx as nx
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("graphs", nargs="+",
                        help="other graph(s) to compare to baseline")
    parser.add_argument("out_dir", type=str,
                        help="directory to which results will be written.")
    parser.add_argument("--kmer", type=int, default=10,
                        help="kmer size for comparison")
    parser.add_argument("--edge_max", type=int, default=0,
                        help="edge-max parameter for vg kmer index")    
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite existing files")                        
    args = args[1:]

    return parser.parse_args(args)

def index_path(graph, options):
    """ get the path of the index given the graph
    """
    return graph + ".index"

def compute_kmer_index(job, graph, options):
    """ run vg index (if necessary) and vg compare on the input
    vg indexes are just created in place, ie same dir as graph,
    so need to have write permission there
    """
    out_index_path = index_path(graph, options)
    do_index = options.overwrite or not os.path.exists(out_index_path)

    index_opts = "-k {}".format(options.kmer)
    if options.edge_max > 0:
        index_opts += " -e {}".format(options.edge_max)

    if do_index:
        os.system("vg index {} {}".format(index_opts, graph))

def comp_path(graph1, graph2, options):
    """ get the path for json output of vg compare
    """
    return os.path.join(options.out_dir, "compare",
                        os.path.splitext(os.path.basename(graph1))[0] + "_vs_" +
                        os.path.splitext(os.path.basename(graph2))[0] + ".json")

def mat_path(options):
    """ get the path of the distance matrix
    """
    return os.path.join(options.out_dir, "distmat.tsv")

def tree_path(options):
    """ path for newick tree
    """
    return os.path.join(options.out_dir, "tree.newick")

def draw_len(weight):
    """ actual weights are between 0 and 1 but vary by many orders of
    magnitude.  try to map them into something for graphviz edge length hint
    """
    if weight < 0.0001:
        return .5
    elif weight < 0.001:
        return .75
    elif weight < 0.01:
        return 1.
    elif weight < 0.1:
        return 1.25
    elif weight < 0.2:
        return 1.6
    else:
        return 1.6 + weight

def cluster_comparisons(options):
    """ scape the comparison files into a distance matrix, cluster into
    a tree
    """    
    mat = dict()
    for graph in options.graphs:
        mat[graph] = dict()

    # make distance matrix in memory
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if graph1 <= graph2:
                jpath = comp_path(graph1, graph2, options)
                jaccard = -1.
                with open(jpath) as f:
                    j = json.loads(f.read())
                    if float(j["union"]) == 0:
                        jaccard = 2.
                    else:
                        jaccard = float(j["intersection"]) / float(j["union"])
                mat[graph1][graph2] = 1. - jaccard
                mat[graph2][graph1] = 1. - jaccard

    # save matrix to file
    with open(mat_path(options), "w") as mat_file:
        mat_buf = "\t" + "\t".join([os.path.basename(i) for i in options.graphs]) + "\n"
        for graph1 in options.graphs:
            mat_buf += os.path.basename(graph1) + "\t"
            mat_buf += "\t".join([str(mat[graph1][i]) for i in options.graphs]) + "\n"
        robust_makedirs(os.path.dirname(mat_path(options)))
        mat_file.write(mat_buf)

    # oops, convert to biopython matrix
    names = [os.path.splitext(os.path.basename(i))[0] for i in options.graphs]
    matrix = []
    for i in xrange(len(options.graphs)):
        row = []
        for j in xrange(i + 1):
            # tree constructor writes 0-distances as 1s for some reason
            # so we hack around here
            val = float(mat[options.graphs[i]][options.graphs[j]])
            if val == 0.:
                val = 1e-10
            elif val == 1.:
                val = 1.1
            row.append(val)
        matrix.append(row)
    dm = _DistanceMatrix(names, matrix)

    # upgma tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    robust_makedirs(os.path.dirname(tree_path(options)))
    Phylo.write(tree, tree_path(options), "newick")

    # png tree -- note : doesn't work in toil
    def f(x):
        if "Inner" in str(x):
            return ""
        else:
            return x
    Phylo.draw_graphviz(tree, label_func = f, node_size=1000, node_shape="s", font_size=10)
    pylab.savefig(tree_path(options).replace("newick", "png"))

    # graphviz
    # get networkx graph
    nxgraph = Phylo.to_networkx(tree)
    # make undirected
    nxgraph = nx.Graph(nxgraph)
    # push names to name labels
    nxgraph = nx.convert_node_labels_to_integers(nxgraph, label_attribute="label")
    for node_id in nxgraph.nodes():
        node = nxgraph.node[node_id]
        if "Inner" in str(node["label"]):
            node["label"] = "\"\""
            node["width"] = 0.001
            node["height"] = 0.001
        else:
            node["fontsize"] = 18
    for edge_id in nxgraph.edges():
        edge = nxgraph.edge[edge_id[0]][edge_id[1]]
        # in graphviz, weight means something else, so make it a label
        weight = float(edge["weight"])
        # undo hack from above
        if weight > 1:
            weight = 1.
        if weight <= 1e-10 or weight == 1.:
            weight = 0.
        edge["weight"] = None
        edge["label"] = "{0:.3g}".format(float(weight) * 100.)
        edge["fontsize"] = 14
        edge["len"] = draw_len(weight)
    nx.write_dot(nxgraph, tree_path(options).replace("newick", "dot"))
    

def compute_kmer_comparison(job, graph1, graph2, options):
    """ run vg compare between two graphs
    """
    out_path = comp_path(graph1, graph2, options)
    graph1_index_path = index_path(graph1, options)
    assert os.path.exists(graph1_index_path)
    graph2_index_path = index_path(graph2, options)
    assert os.path.exists(graph2_index_path)

    do_comp = options.overwrite or not os.path.exists(out_path)
    
    if do_comp:
        robust_makedirs(os.path.dirname(out_path))        
        os.system("vg compare {} {} > {}".format(graph1, graph2, out_path))

def compute_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    for graph1 in options.graphs:
        for graph2 in options.graphs:
            if graph1 <= graph2:
                job.addChildJobFn(compute_kmer_comparison, graph1, graph2, options)

def compute_kmer_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """
    # do all the indexes
    for graph in options.graphs:
        job.addChildJobFn(compute_kmer_index, graph, options)

    # do the comparisons
    job.addFollowOnJobFn(compute_comparisons, options)
    
def main(args):
    
    options = parse_args(args) 
    
    RealTimeLogger.start_master()

    for graph in options.graphs:
        if os.path.splitext(graph)[1] != ".vg":
            raise RuntimeError("Input graphs expected to have .vg extension")

    # Make a root job
    root_job = Job.wrapJobFn(compute_kmer_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()

    # Do the drawing outside toil to get around weird import problems
    cluster_comparisons(options)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

