#!/usr/bin/env python2.7
"""
parallelMappingEvaluation.py: Run the mapping evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading

from toil.job import Job
import toillib
from toillib import RealTimeLogger

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("server_list", type=argparse.FileType("r"),
        help="TSV file continaing <region>\t<url> lines for servers to test")
    parser.add_argument("sample_store", type=toillib.IOStore.get,
        help="sample input IOStore with <region>/<sample>/<sample>.bam.fq")
    parser.add_argument("out_store", type=toillib.IOStore.get,
        help="output IOStore to create and fill with alignments and stats")
    parser.add_argument("--server_version", default="v0.6.g",
        help="server version to add to URLs")
    parser.add_argument("--sample_limit", type=int, default=1, 
        help="number of samples to use")
    parser.add_argument("--edge_max", type=int, default=0, 
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=10, 
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--bin_url",
        default="https://hgwdev.sdsc.edu/~anovak/bin",
        help="URL to download sg2vg and vg binaries from, without Docker")
    parser.add_argument("--ca_cert_url",
        default="https://www.incommon.org/certificates/repository/"
        "sha384%20Intermediate%20cert.txt",
        help="CA PEM certificate to download and use to get the binaries")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    

def run_all_alignments(job, options):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    
    # Get the CA file needed to verify the binary downloads
    # The certificate has to come from somewhere trusted by system wget.
    RealTimeLogger.get().info("Retrieving CA certificate {}".format(
        options.ca_cert_url))
    ca_file = "{}/cert.pem".format(job.fileStore.getLocalTempDir())
    subprocess.check_call(["wget", options.ca_cert_url,
        "-O", ca_file])
    
    # Fix up any \r line endings
    RealTimeLogger.get().info("Rewriting certificate")
    subprocess.check_call(["tr", "\\r", "\\n"], stdin=open(ca_file, "r"),
        stdout=open(ca_file + ".fixed", "w"))
    os.rename(ca_file + ".fixed", ca_file)
    
    
    # Retrieve binaries we need, using the certificate
    RealTimeLogger.get().info("Retrieving binaries from {}".format(
        options.bin_url))
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    toillib.robust_makedirs(bin_dir)
    subprocess.check_call(["wget", "{}/sg2vg".format(options.bin_url),
        "--ca-certificate", ca_file, "-O", "{}/sg2vg".format(bin_dir)])
    subprocess.check_call(["wget", "{}/vg".format(options.bin_url),
        "--ca-certificate", ca_file, "-O", "{}/vg".format(bin_dir)])
        
    # Make them executable
    os.chmod("{}/sg2vg".format(bin_dir), 0o744)
    os.chmod("{}/vg".format(bin_dir), 0o744)
    
    # Upload the bin directory to the file store
    bin_dir_id = toillib.write_global_directory(job.fileStore, bin_dir,
        cleanup=True)
    
    # Make sure we skip the header
    is_first = True
    
    for line in options.server_list:
        if is_first:
            # This is the header, skip it.
            is_first = False
            continue
        
        # We need to read each non-header line
        
        # Break it into its fields
        parts = line.split("\t")
        
        if parts[0].startswith("#"):
            # Skip comments
            continue
            
        if parts[0].startswith("\n"):
            # Skip newlines
            continue
            
        # Pull out the first 3 fields
        region, url, generator = parts[0:3]
        
        # We cleverly just split the lines out to different nodes
        job.addChildJobFn(run_region_alignments, options, bin_dir_id, region,
            url, cores=16, memory="100G", disk="50G")
            
        # Say what we did
        RealTimeLogger.get().info("Running child for {}".format(parts[1]))
        

def run_region_alignments(job, options, bin_dir_id, region, url):
    """
    For the given region, download, index, and then align to the given graph.
    
    """
    
    RealTimeLogger.get().info("Running on {} for {}".format(url, region))
    
    # Download the binaries
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    toillib.read_global_directory(job.fileStore, bin_dir_id, bin_dir)
    
    
    # Get graph basename (last URL component) from URL
    basename = re.match(".*/(.*)/$", url).group(1)
        
    # Get graph name (without region and its associated dash) from basename
    graph_name = basename.replace("-{}".format(region), "").replace(
        "{}-".format(region), "")
    
    # Make the real URL with the version
    versioned_url = url + options.server_version
    
    # Work out where the graph goes
    # it will be graph.vg in here
    graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
    toillib.robust_makedirs(graph_dir)
    
    graph_filename = "{}/graph.vg".format(graph_dir)
    
    # Download and fix up the graph with this ugly subprocess pipeline
    # sg2vg "${URL}" -u | vg view -Jv - | vg mod -X 100 - | vg ids -s - > "graphs/${BASENAME}.vg"
    
    with open(graph_filename, "w") as output_file:
    
        RealTimeLogger.get().info("Downloading {} to {}".format(versioned_url,
            graph_filename))
    
        # Hold all the popen objects we need for this
        tasks = []
        
        # Do the download
        tasks.append(subprocess.Popen(["{}/sg2vg".format(bin_dir),
            versioned_url, "-u"], stdout=subprocess.PIPE))
        
        # Pipe through zcat
        tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "view", "-Jv",
            "-"], stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
        
        # And cut
        tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "mod", "-X100",
            "-"], stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
            
        # And uniq
        tasks.append(subprocess.Popen(["{}/vg".format(bin_dir), "ids", "-s",
            "-"], stdin=tasks[-1].stdout, stdout=output_file))
            
        # Did we make it through all the tasks OK?
        for task in tasks:
            if task.wait() != 0:
                raise RuntimeError("Pipeline step returned {}".format(
                    task.returncode))
                    
    # Now run the indexer.
    # TODO: support both indexing modes
    RealTimeLogger.get().info("Indexing {}".format(graph_filename))
    subprocess.check_call(["{}/vg".format(bin_dir), "index", "-s", "-k",
        str(options.kmer_size), "-e", str(options.edge_max), graph_filename])
        
    # Now save the indexed graph directory to the file store. It can be cleaned
    # up since only our children use it.
    index_dir_id = toillib.write_global_directory(job.fileStore, graph_dir,
        cleanup=True)
                    
    # Where do we look for samples for this region in the input?
    region_dir = region.upper()
    
    # What samples do we do? List input sample names up to the given limit.
    input_samples = list(options.sample_store.list_input_directory(
        region_dir))[:options.sample_limit]
    
    # Work out the directory for the alignments to be dumped in in the output
    alignment_dir = "alignments/{}/{}".format(region, graph_name)
    
    # Also for statistics
    stats_dir = "stats/{}/{}".format(region, graph_name)
    
    
    for sample in input_samples:
        # Split out over each sample
        
        RealTimeLogger.get().info("Queueing alignment of {} to {} {}".format(
            sample, graph_name, region))
    
        # For each sample, know the FQ name
        sample_fastq = "{}/{}/{}.bam.fq".format(region_dir, sample, sample)
        
        # And know where we're going to put the output
        alignment_file_key = "{}/{}.gam".format(alignment_dir, sample)
        stats_file_key = "{}/{}.json".format(stats_dir, sample)
    
        # Go and bang that input fastq against the correct indexed graph.
        # Its output will go to the right place in the output store.
        job.addChildJobFn(run_alignment, options, bin_dir_id, region,
            index_dir_id, sample_fastq, alignment_file_key, stats_file_key,
            cores=16, memory="100G", disk="50G")
   
def run_alignment(job, options, bin_dir_id, region, index_dir_id,
    sample_fastq_key, alignment_file_key, stats_file_key):
    """
    Align the the given fastq from the input store against the given indexed
    graph (in the file store as a directory) and put the GAM and statistics in
    the given output keys in the output store.
    
    """
    
    # Download the binaries
    bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
    toillib.read_global_directory(job.fileStore, bin_dir_id, bin_dir)
    
    # Download the indexed graph to a directory we can use
    graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
    toillib.read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)
    
    # Also we need the sample fastq
    fastq_file = "{}/input.fq".format(job.fileStore.getLocalTempDir())
    options.sample_store.read_input_file(sample_fastq_key, fastq_file)
    
    # And temp files for our aligner output and stats
    output_file = "{}/output.gam".format(job.fileStore.getLocalTempDir())
    stats_file = "{}/stats.json".format(job.fileStore.getLocalTempDir())
    
    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:
    
        # Start the aligner and have it write to the file
        
        # Plan out what to run
        vg_parts = ["{}/vg".format(bin_dir), "map", "-f", fastq_file, "-i",
            "-n3", "-M2", "-k", str(options.kmer_size), graph_file]
        
        RealTimeLogger.get().info("Running VG: {}".format(" ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        process = subprocess.Popen(vg_parts, stdout=alignment_file)
            
        if process.wait() != 0:
            # Complain if vg dies
            raise RuntimeError("vg died with error {}".format(
                process.returncode))
                
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        
                
    RealTimeLogger.get().info("Aligned {}".format(output_file))
           
    # Read the alignments in in JSON-line format
    view = subprocess.Popen(["{}/vg".format(bin_dir), "view", "-aj",
        output_file], stdout=subprocess.PIPE)
       
    # Count up the stats: total reads, total mapped at all, total multimapped,
    # primary alignment score counts, secondary alignment score counts, and
    # aligner run time in seconds.
    
    stats = {
        "total_reads": 0,
        "total_mapped": 0,
        "total_multimapped": 0,
        "primary_scores": collections.Counter(),
        "primary_mismatches": collections.Counter(),
        "secondary_scores": collections.Counter(),
        "secondary_mismatches": collections.Counter(),
        "run_time": run_time
    }
        
    for line in view.stdout:
        # Parse the alignment JSON
        alignment = json.loads(line)
        
        if alignment.has_key("score"):
            # This alignment is aligned.
            # Grab its score
            score = alignment["score"]
        
            # Calculate the mismatches
            length = len(alignment["sequence"])
            matches = 0
            for mapping in alignment.get("path", {}).get("mapping", []):
                for edit in mapping.get("edit", []):
                    if (not edit.has_key("sequence") and 
                        edit.get("to_length", None) == edit.get(
                        "from_length", None)):
                        
                        # We found a perfect match edit. Grab its length
                        matches += edit["from_length"]
                        
            # Calculate mismatches as what's left
            mismatches = length - matches
                    
        
            if alignment.get("is_secondary", False):
                # It's a multimapping. We can have max 1 per read, so it's a
                # multimapped read.
                
                # Log its stats as multimapped
                stats["total_multimapped"] += 1
                stats["secondary_scores"][score] += 1
                stats["secondary_mismatches"][mismatches] += 1
            else:
                # Log its stats as primary. We'll get exactly one of these per
                # read with any mappings.
                stats["total_mapped"] += 1
                stats["primary_scores"][score] += 1
                stats["primary_mismatches"][mismatches] += 1
                
                # We won't see an unaligned primary alignment for this read, so
                # count the read
                stats["total_reads"] += 1
        
        elif not alignment.get("is_secondary", False):
            # We have an unmapped primary "alignment"
            
            # Count the read by its primary alignment
            stats["total_reads"] += 1
                
    with open(stats_file, "w") as stats_handle:
        # Save the stats as JSON
        json.dump(stats, stats_handle)
        
    # Now send the output files (alignment and stats) to the output store where
    # they belong.
    options.out_store.write_output_file(output_file, alignment_file_key)
    options.out_store.write_output_file(stats_file, stats_file_key)
    
        
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    RealTimeLogger.start_master()
    
    # Pre-read the input file so we don't try to send file handles over the
    # network.
    options.server_list = list(options.server_list)
    
    # Make a root job
    root_job = Job.wrapJobFn(run_all_alignments, options,
        cores=16, memory="100G", disk="50G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

