#!/usr/bin/env python2.7
"""
Run vg compare on a bunch of graphs and make a table on how there kmer sets match up. 
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from toil.job import Job
from parallelMappingEvaluation import RealTimeLogger, robust_makedirs
from callVariants import call_path, augmented_vg_path, alignment_read_tag, alignment_map_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path 

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("baseline", type=str,
                        help="baseline graph to compare all other graphs to")
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files. (must have been run through callVariants.py!)")
    parser.add_argument("--out_dir", type=str, default="variants",
                        help="output directory (used for callVariants.py)")
    parser.add_argument("--graph_dir", type=str, default="graphs",
                        help="name of input graphs directory")
    parser.add_argument("--index_ext", type=str, default=".index",
                        help="extension to find input grpah index")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite files if dont exist")
    parser.add_argument("out_tsv", type=str,
                        help="path for results table file (tab-separated)")
    parser.add_argument("--kmer", type=int, default=10,
                        help="kmer size for comparison")
    parser.add_argument("--edge_max", type=int, default=0,
                        help="edge-max parameter for vg kmer index")
                        
    args = args[1:]
        
    return parser.parse_args(args)

def index_path(graph, options):
    """ get the path of the index given the graph
    """
    return graph + ".index"

def comparison_path(baseline, graph, options):
    """ get a unique path, in the temporary folder, for
    json comparison output of two graphs.  note that graph
    filename may not be unique, so we extract some tags out of its path
    """
    bname = os.path.splitext(os.path.basename(baseline))[0]
    gname = os.path.splitext(os.path.basename(graph))[0]
    try:
        read_tag = alignment_read_tag(graph, options)
    except:
        read_tag = ""
    try:
        map_tag = alignment_map_tag(graph, options)
    except:
        map_tag = ""

    tempdir = os.path.join(options.out_dir, "temp")
    
    return os.path.join(tempdir, "{}{}{}{}.json".format(
        bname, map_tag, read_tag, gname))
    
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
        RealTimeLogger.get().info("graph path {}".format(graph))
        assert os.path.isfile(graph)
        os.system("vg index {} {}".format(index_opts, graph))

def compute_comparison(job, baseline, graph, options):
    """ run vg compare between two graphs
    """
    graph_index_path = index_path(graph, options)
    assert os.path.exists(graph_index_path)
    baseline_index_path = index_path(baseline, options)
    assert os.path.exists(baseline_index_path)

    out_path = comparison_path(baseline, graph, options)
    do_comp = options.overwrite or not os.path.exists(out_path)
    
    if do_comp:        
        os.system("vg compare {} {} > {}".format(baseline, graph, out_path))
           
def count_gam_paths(graph, options):
    """ get number of snps by counting the paths in an alignment. only
    workds for agumented vg graphs made by the callVariants script
    """
    aug_tag = os.path.basename(augmented_vg_path("/a/b/c", options))[1:]
    call_tag = os.path.basename(call_path("/a/b/c", options))[1:]
    if graph.find(aug_tag) < 0:
        return -1
    gam = graph.replace(aug_tag, call_tag)
    cmd = "vg view -a -j {} | jq .path | jq length | wc -l".format(gam)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def jaccard_dist(comp_json):
    """ get a distance from the kmer set comparison json output
    """
    intersection_size = float(comp_json["intersection"])
    union_size = float(comp_json["union"])
    if union_size == 0.:
        return -1
    else:
        return 1. - intersection_size / union_size
    
def summarize_comparisons(job, options):
    """ make the table by scraping together all the comparison
    json files
    """
    table =  "#\t{}\t\t\t\t\n".format(os.path.splitext(os.path.basename(options.baseline))[0])
    table += "#graph\tgraph_dist\tlinear_dist\taugmented_dist\t\delta_linear\t\delta_augmented\n"
    for gam in options.in_gams:
        comp_path = comparison_path(options.baseline, graph_path(gam, options), options)
        with open(comp_path) as f:
            j = json.loads(f.read())
            graph_dist = jaccard_dist(j)
        aug_comp_path = comparison_path(options.baseline, augmented_vg_path(gam, options), options)
        with open(aug_comp_path) as f:
            j = json.loads(f.read())
            aug_graph_dist = jaccard_dist(j)
        lin_comp_path = comparison_path(options.baseline, linear_vg_path(gam, options), options)
        with open(lin_comp_path) as f:
            j = json.loads(f.read())
            lin_graph_dist = jaccard_dist(j)

        delta_lin = lin_graph_dist - graph_dist
        delta_aug = aug_graph_dist - graph_dist

        table += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
            os.path.splitext(os.path.basename(graph_path(gam, options)))[0],
            graph_dist,
            lin_graph_dist,
            aug_graph_dist,
            delta_lin,
            delta_aug)

    with open(options.out_tsv, "w") as ofile:
        ofile.write(table)

def compute_all_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    for gam in options.in_gams:
        job.addChildJobFn(compute_comparison, options.baseline,
                          graph_path(gam, options), options)
        job.addChildJobFn(compute_comparison, options.baseline,
                          augmented_vg_path(gam, options), options)
        job.addChildJobFn(compute_comparison, options.baseline,
                          linear_vg_path(gam, options), options)

    job.addFollowOnJobFn(summarize_comparisons, options)
                
def compute_all_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """

    # do all the indexes
    job.addChildJobFn(compute_kmer_index, options.baseline, options, cores=1)
    for gam in options.in_gams:
        if graph_path(gam, options) != options.baseline:
            job.addChildJobFn(compute_kmer_index, graph_path(gam, options), options)
        if augmented_vg_path(gam, options) != options.baseline:
            job.addChildJobFn(compute_kmer_index, augmented_vg_path(gam, options), options)
        if linear_vg_path(gam, options) != options.baseline:
            job.addChildJobFn(compute_kmer_index, linear_vg_path(gam, options), options)

    # do the comparisons
    job.addFollowOnJobFn(compute_all_comparisons, options)

    
def main(args):
    
    options = parse_args(args) 
    
    RealTimeLogger.start_master()

    for gam in options.in_gams:
        if len(gam.split("/")) < 3 or os.path.splitext(gam)[1] != ".gam":
            raise RuntimeError("Input gam paths must be of the form "
                               ".../<alg>/<reads>/<filename>.gam")
    # Make a root job
    root_job = Job.wrapJobFn(compute_all_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

