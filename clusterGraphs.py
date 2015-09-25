#!/usr/bin/env python2.7
"""
Cluster some graphs (using vg compare for pairwise distances) based on their similarity.
Current implementation : neighbour joining tree using Jaccard distance matrix
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
from toil.job import Job
from parallelMappingEvaluation import RealTimeLogger, robust_makedirs

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

    if do_index:
        os.system("vg index -k {} {}".format(options.kmer, graph))

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

def cluster_comparisons(job, options):
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
            row.append(mat[options.graphs[i]][options.graphs[j]])
        matrix.append(row)
    dm = _DistanceMatrix(names, matrix)

    # upgma tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    robust_makedirs(os.path.dirname(tree_path(options)))
    Phylo.write(tree, tree_path(options), "newick")

    # png tree
    try:
        import pylab
        Phylo.draw_graphviz(tree)
        pylab.savefig(tree_path(options).replace("newick", "png"))
    except:
        pass

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

    job.addFollowOnJobFn(cluster_comparisons, options)
                
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
    
    RealTimeLogger.start_master(options)

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
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

