#!/usr/bin/env python2.7
"""
Run vg compare on a bunch of graphs and make a table on how there kmer sets match up. 
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from operator import sub
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs
from callVariants import call_path, augmented_vg_path, alignment_read_tag, alignment_map_tag
from callVariants import graph_path, index_path, augmented_vg_path, linear_vg_path, linear_vcf_path
from callVariants import alignment_region_tag

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("baseline", type=str,
                        help="baseline graph to compare all other graphs to. "
                        "(just the name part, ex snp1000g. graph_dir/<baseline>-<region>.vg"
                        " must exist for each region)")
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

def compare_out_path(options):
    """ get root output dir for comparison output
    """
    return os.path.join(options.out_dir, "compare")

def json_out_path(options):
    """ where to put the json compareison output
    """
    return os.path.join(compare_out_path(options), "json")

def dist_tsv_path(options):
    """ path where the tsv for tdistance goes
    """
    return os.path.join(compare_out_path(options),
                        "call_dist_{}.tsv".format(options.baseline))

def acc_tsv_path(options):
    """ path where the tsv for prec/rec/f1 goes
    """
    return os.path.join(compare_out_path(options),
                        "call_acc_{}.tsv".format(options.baseline))

def count_tsv_path(options):
    return os.path.join(compare_out_path(options),
                        "call_count_{}.tsv".format(options.baseline))

def baseline_path(gam, options):
    """ put together path for baseline graph 
    """
    region = alignment_region_tag(gam, options)
    return os.path.join(options.graph_dir,
                        "{}-{}.vg".format(options.baseline, region))

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
    
    return os.path.join(json_out_path(options), "{}{}{}{}.json".format(
        bname, map_tag, read_tag, gname))
    
def compute_kmer_index(job, graph, options):
    """ run vg index (if necessary) and vg compare on the input
    vg indexes are just created in place, ie same dir as graph,
    so need to have write permission there
    """
    out_index_path = index_path(graph, options)
    do_index = options.overwrite or not os.path.exists(out_index_path)

    index_opts = "-s -k {}".format(options.kmer)
    if options.edge_max > 0:
        index_opts += " -e {}".format(options.edge_max)
    
    if do_index:
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
           
def count_gam_paths(gam, options):
    """ get number of snps by counting the paths in an alignment. only
    workds for agumented vg graphs made by the callVariants script
    """
    if not os.path.exists(gam):
        return -1
    cmd = "vg view -a -j {} | jq .path | jq length | wc -l".format(gam)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def count_vcf_snps(vcf, options):
    """ get number of snps from bcftools
    """
    if not os.path.exists(vcf):
        return -1
    cmd = "./vcfCountSnps.sh {}".format(vcf)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    num_paths = int(output.strip())
    return num_paths

def jaccard_dist(comp_json_path):
    """ get a distance from the kmer set comparison json output
    """
    if not os.path.isfile(comp_json_path):
        return -1.
    with open(comp_json_path) as f:
        comp_json = json.loads(f.read())

    intersection_size = float(comp_json["intersection"])
    union_size = float(comp_json["union"])
    if union_size == 0.:
        return -1.
    else:
        return 1. - intersection_size / union_size

def accuracy(comp_json_path):
    """ compute precision, recall, f1 from the kmer set comparison json output
    (assuming db1 is the "truth")
    """
    if not os.path.isfile(comp_json_path):
        return -1., -1., -1.
    with open(comp_json_path) as f:
        comp_json = json.loads(f.read())
        
    intersection_size = float(comp_json["intersection"])
    db1_size = float(comp_json["db1_total"])
    db2_size = float(comp_json["db2_total"])
    if db2_size > 0:
        precision = intersection_size / db2_size
    else:
        precision = -1.
    if db1_size > 0:
        recall = intersection_size / db1_size
    else:
        recall = -1.
    if recall >= 0 and precision >= 0 and precision + recall > 0:
        f1 = 2. * ((precision * recall) / (precision + recall))
    else:
        f1 = -1.
    return precision, recall, f1

def dist_table(options):
    """ make the jaccard distance table by scraping together all the comparison
    json files
    """
    # tsv header
    dist_table =  "#\t{}\t\t\t\t\n".format(options.baseline)
    dist_table += "#graph\tgraph_dist\tlinear_dist\taugmented_dist\tdelta_linear\tdelta_augmented\n"
    
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        comp_path = comparison_path(baseline, graph_path(gam, options), options)
        graph_dist = jaccard_dist(comp_path)
        aug_comp_path = comparison_path(baseline, augmented_vg_path(gam, options), options)
        aug_graph_dist = jaccard_dist(aug_comp_path)
        lin_comp_path = comparison_path(baseline, linear_vg_path(gam, options), options)
        lin_graph_dist = jaccard_dist(lin_comp_path)

        delta_lin = lin_graph_dist - graph_dist
        delta_aug = aug_graph_dist - graph_dist

        dist_table += "{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n".format(
            os.path.splitext(os.path.basename(graph_path(gam, options)))[0],
            graph_dist,
            lin_graph_dist,
            aug_graph_dist,
            delta_lin,
            delta_aug)

    with open(dist_tsv_path(options), "w") as ofile:
        ofile.write(dist_table)


def acc_table(options):
    """ make the accuracy table by scraping together all the comparison
    json files
    """
    # tsv header
    acc_table =  "#\t{}\t\t\t\t\t\t\t\t\t\t\t\t\n".format(options.baseline)
    acc_table += "#graph\tgraph_prec\tgraph_rec\tgraph_f1"
    acc_table += "\tlinear_prec\tlinear_rec\tlinear_f1"
    acc_table += "\taugmented_prec\taugmented_rec\taugmented_f1"
    acc_table += "\tprec_delta_linear\trec_delta_linear\tf1_delta_linear"
    acc_table += "\tprec_delta_augmented\trec_delta_augmented\tf1_delta_augmented\n"
    
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        comp_path = comparison_path(baseline, graph_path(gam, options), options)
        graph_acc = accuracy(comp_path)
        aug_comp_path = comparison_path(baseline, augmented_vg_path(gam, options), options)
        aug_graph_acc = accuracy(aug_comp_path)
        lin_comp_path = comparison_path(baseline, linear_vg_path(gam, options), options)
        lin_graph_acc = accuracy(lin_comp_path)

        delta_lin_acc = map(sub, lin_graph_acc, graph_acc)
        delta_aug_acc = map(sub, aug_graph_acc, graph_acc)

        acc_table += "{}\t{:.4}\t{:.4}\t{:.4}\t".format(
            os.path.splitext(os.path.basename(graph_path(gam, options)))[0],
            graph_acc[0], graph_acc[1], graph_acc[2])
        acc_table += "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t".format(
            lin_graph_acc[0], lin_graph_acc[1], lin_graph_acc[2],
            aug_graph_acc[0], aug_graph_acc[1], aug_graph_acc[2])
        acc_table +="{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n".format(
            delta_lin_acc[0], delta_lin_acc[1], delta_lin_acc[2],
            delta_aug_acc[0], delta_aug_acc[1], delta_aug_acc[2])

    with open(acc_tsv_path(options), "w") as ofile:
        ofile.write(acc_table)

def snp_count_table(options):
    """ make a table of snp counts.  there are serious problems with this now:
    1) don't have snp count for baseline (as it's not gam or vcf)
    2) snps counted differenty for gam/vcf (multiple alternates at same site
    counted in former but not latter)
    """
    # tsv header
    count_table =  "#\t{}\t\n".format(options.baseline)
    count_table += "#graph\tlinear_snp_count\taugmented_snp_count\n"

    for gam in options.in_gams:
        linear_vcf = linear_vcf_path(gam, options) + ".gz"
        vg_gam = call_path(gam, options)
        vcf_snps = count_vcf_snps(linear_vcf, options)
        gam_snps = count_gam_paths(vg_gam, options)

        count_table +="{}\t{}\t{}\n".format(
            graph_path(gam, options),
            vcf_snps,
            gam_snps)

    with open(count_tsv_path(options), "w") as ofile:
        ofile.write(count_table)
    
def compute_all_comparisons(job, options):
    """ run vg compare in parallel on all the graphs,
    outputting a json file for each
    """
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        job.addChildJobFn(compute_comparison, baseline,
                          graph_path(gam, options), options)
        job.addChildJobFn(compute_comparison, baseline,
                          augmented_vg_path(gam, options), options)
        job.addChildJobFn(compute_comparison, baseline,
                          linear_vg_path(gam, options), options)

def compute_all_indexes(job, options):
    """ run everything (root toil job)
    first all indexes are computed,
    then all comparisons (follow on)
    then summary (follow on of that)
    """

    # do all the indexes
    baseline_set = set()
        
    for gam in options.in_gams:
        baseline = baseline_path(gam, options)
        if not os.path.isfile(baseline):
            raise RuntimeError("baseline {} for gam {} not found".format(baseline, gam))
        if baseline not in baseline_set:
            job.addChildJobFn(compute_kmer_index, baseline, options, cores=1)
            baseline_set.add(baseline)
        if graph_path(gam, options) != baseline:
            job.addChildJobFn(compute_kmer_index, graph_path(gam, options), options)
        if augmented_vg_path(gam, options) != baseline:
            job.addChildJobFn(compute_kmer_index, augmented_vg_path(gam, options), options)
        if linear_vg_path(gam, options) != baseline:
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
    robust_makedirs(json_out_path(options))
    robust_makedirs(compare_out_path(options))
                    
    # Make a root job
    root_job = Job.wrapJobFn(compute_all_indexes, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()

    # make some tables from the json comparison output
    dist_table(options)
    acc_table(options)
    snp_count_table(options)
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

