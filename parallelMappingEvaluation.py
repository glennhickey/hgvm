#!/usr/bin/env python2.7
"""
parallelMappingEvaluation.py: Run the mapping evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections
import logging, logging.handlers, SocketServer, struct, socket, threading

from toil.job import Job

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
    parser.add_argument("sample_dir",
        help="sample input directory with <region>/<sample>/<sample>.bam.fq")
    parser.add_argument("out_dir",
        help="output directory to create and fill with alignments and stats")
    parser.add_argument("--server_version", default="v0.6.g",
        help="server version to add to URLs")
    parser.add_argument("--sample_limit", type=int, default=1, 
        help="number of samples to use")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
class LoggingDatagramHandler(SocketServer.DatagramRequestHandler):
    """
    Receive logging messages from the jobs and display them on the master.
    
    Uses length-prefixed JSON message encoding.
    """
    
    def handle(self):
        """
        Handle messages coming in over self.connection.
        
        Messages are 4-byte-length-prefixed JSON-encoded logging module records.
        """
        
        while True:
            # Loop until we run out of messages
        
            # Parse the length
            length_data = self.rfile.read(4)
            if len(length_data) < 4:
                # The connection was closed, or we didn't get enough data
                # TODO: complain?
                break
                
            # Actually parse the length
            length = struct.unpack(">L", length_data)[0]
            
            # This is where we'll put the received message
            message_parts = []
            length_received = 0
            while length_received < length:
                # Keep trying to get enough data
                part = self.rfile.read(length - length_received)
                
                length_received += len(part)
                message_parts.append(part)
                
            # Stitch it all together
            message = "".join(message_parts)

            try:
            
                # Parse it as JSON
                message_attrs = json.loads(message)
                
                # Fluff it up into a proper logging record
                record = logging.makeLogRecord(message_attrs)
            except:
                logging.error("Malformed record")
                
            # TODO: do log level filtering
            logging.getLogger("remote").handle(record)
            
class JSONDatagramHandler(logging.handlers.DatagramHandler):
    """
    Send logging records over UDP serialized as JSON.
    """
    
    def makePickle(self, record):
        """
        Actually, encode the record as length-prefixed JSON instead.
        """
        
        json_string = json.dumps(record.__dict__)
        length = struct.pack(">L", len(json_string))
        
        return length + json_string
        
class RealTimeLogger(object):
    """
    All-static class for getting a logger that logs over UDP to the master.
    """
    
    # We keep the master host and port as class data
    master_host = None
    master_port = None
    
    # Also the logger
    logger = None
    
    # The master keeps a server and thread
    logging_server = None
    server_thread = None
  
    @classmethod
    def start_master(cls, options):
        """
        Start up the master server and put its details into the options
        namespace.
        
        """
        
        logging.basicConfig(level=logging.DEBUG)
    
        # Start up the logging server
        cls.logging_server = SocketServer.ThreadingUDPServer(("0.0.0.0", 0),
            LoggingDatagramHandler)
            
        # Set up a thread to do all the serving in the background and exit when we
        # do
        cls.server_thread = threading.Thread(
            target=cls.logging_server.serve_forever)
        cls.server_thread.daemon = True
        cls.server_thread.start()
        
        # Set options for logging in the class and the options namespace
        cls.master_host = socket.getfqdn()
        options.log_host = cls.master_host
        cls.master_port = cls.logging_server.server_address[1]
        options.log_port = cls.master_port
        
    @classmethod
    def stop_master(cls):
        """
        Stop the server on the master.
        
        """
        
        cls.logging_server.shutdown()
        cls.server_thread.join()
  
    @classmethod
    def set_master(cls, options):
        """
        Set the master info.
        
        """
        
        cls.master_host = options.log_host
        cls.master_port = options.log_port
        
        
    @classmethod
    def get(cls):
        """
        Get the logger that logs to master.
        
        Note that if the master logs here, you will see the message twice,
        since it still goes to the normal log handlers too.
        """
        
        if cls.logger is None:
            # Only do the setup once, so we don't add a handler every time we
            # log
            cls.logger = logging.getLogger('realtime')
            cls.logger.setLevel(logging.DEBUG)
            cls.logger.addHandler(JSONDatagramHandler(cls.master_host,
                cls.master_port))
        
        return cls.logger

    
def robust_makedirs(directory):
    """
    Make a directory when other nodes may be trying to do the same on a shared
    filesystem.
    
    """

    if not os.path.exists(directory):
        try:
            # Make it if it doesn't exist
            os.makedirs(directory)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass
            
    # Make sure it exists and is a directory
    assert(os.path.exists(directory) and os.path.isdir(directory))
    
def run_all_alignments(job, options):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    
    RealTimeLogger.set_master(options)

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
        job.addChildJobFn(run_region_alignments, options, region, url,
            cores=16, memory="100G", disk="50G")
            
        # Say what we did
        RealTimeLogger.get().info("Running child for {}".format(parts[1]))
        

def run_region_alignments(job, options, region, url):
    """
    For the given region, download, index, and then align to the given graph.
    
    """
    
    RealTimeLogger.set_master(options)
    
    RealTimeLogger.get().info("Running on {} for {}".format(url, region))
    
    # Get graph basename (last URL component) from URL
    basename = re.match(".*/(.*)/$", url).group(1)
        
    # Get graph name (without region and its associated dash) from basename
    graph_name = basename.replace("-{}".format(region), "").replace(
        "{}-".format(region), "")
    
    # Make the real URL with the version
    versioned_url = url + options.server_version
    
    # Work out where the graph goes
    # it will be <graph_name>.vg in here
    graph_dir = "{}/graphs/{}".format(options.out_dir, region)
    robust_makedirs(graph_dir)
    
    graph_filename = "{}/{}.vg".format(graph_dir, graph_name)
    
    # Download and fix up the graph with this ugly subprocess pipeline
    # sg2vg "${URL}" -u | vg view -Jv - | vg mod -X 100 - | vg ids -s - > "graphs/${BASENAME}.vg"
    
    with open(graph_filename, "w") as output_file:
    
        RealTimeLogger.get().info("Downloading {} to {}".format(versioned_url,
            graph_filename))
    
        # Hold all the popen objects we need for this
        tasks = []
        
        # Do the download
        tasks.append(subprocess.Popen(["sg2vg", versioned_url, "-u"],
            stdout=subprocess.PIPE))
        
        # Pipe through zcat
        tasks.append(subprocess.Popen(["vg", "view", "-Jv", "-"],
            stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
        
        # And cut
        tasks.append(subprocess.Popen(["vg", "mod", "-X100", "-"],
            stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
            
        # And uniq
        tasks.append(subprocess.Popen(["vg", "ids", "-s", "-"],
            stdin=tasks[-1].stdout, stdout=output_file))
            
        # Did we make it through all the tasks OK?
        for task in tasks:
            if task.wait() != 0:
                raise RuntimeError("Pipeline step returned {}".format(
                    task.returncode))
                    
    # Now run the indexer.
    # TODO: support both indexing modes
    RealTimeLogger.get().info("Indexing {}".format(graph_filename))
    subprocess.check_call(["vg", "index", "-sk10", graph_filename])
                    
    # Where do we look for samples for this region?
    region_dir = "{}/{}".format(options.sample_dir, region.upper()) 
    
    # Work out the directory for the alignments to be dumped in
    alignment_dir = "{}/alignments/{}/{}".format(options.out_dir, region,
        graph_name)
    robust_makedirs(alignment_dir)
    
    # Also for statistics
    stats_dir = "{}/stats/{}/{}".format(options.out_dir, region,
        graph_name)
    robust_makedirs(stats_dir)
    
    # Split out over each sample
    for sample in list(os.listdir(region_dir))[:options.sample_limit]:
        # For each sample up to the limit
        
        RealTimeLogger.get().info("Queueing alignment of {} to {}".format(
            sample, region))
    
        # For each sample, know the FQ name
        sample_fastq = "{}/{}/{}.bam.fq".format(region_dir, sample, sample)
        
        # And know where we're going to put the output
        alignment_file = "{}/{}.gam".format(alignment_dir, sample)
        stats_file = "{}/{}.json".format(stats_dir, sample)
    
        # Go and bang that fastq against the correct graph.
        # Its output will go to the right place in the output directory.
        job.addChildJobFn(run_alignment, options, region, graph_filename,
            sample_fastq, alignment_file, stats_file, cores=32, memory="240G",
            disk="100G")
   
def run_alignment(job, options, region, graph_file, sample_fastq, output_file,
    stats_file):
    """
    Align the the given fastq against the given graph and put the GAM and
    statistics in the given files.
    
    """
    
    RealTimeLogger.set_master(options)
    
    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:
    
        # Start the aligner and have it write to the file
        
        # Plan out what to run
        vg_parts = ["vg", "map", "-f", sample_fastq, "-i", "-n3", "-M2", "-k10",
            graph_file]
            
        RealTimeLogger.get().info("Running VG: {}".format(" ".join(vg_parts)))
        
        process = subprocess.Popen(vg_parts, stdout=alignment_file)
            
        if process.wait() != 0:
            # Complain if vg dies
            raise RuntimeError("vg died with error {}".format(
                process.returncode))
                
    RealTimeLogger.get().info("Aligned {}".format(output_file))
           
    # Read the alignments in in JSON-line format
    view = subprocess.Popen(["vg", "view", "-aj", output_file],
        stdout=subprocess.PIPE)
       
    # Count up the stats: total reads, total mapped at all, total multimapped,
    # primary alignment score counts, secondary alignment score counts
    
    stats = {
        "total_reads": 0,
        "total_mapped": 0,
        "total_multimapped": 0,
        "primary_perfect": 0,
        "primary_scores": collections.Counter(),
        "primary_mismatches": collections.Counter(),
        "secondary_scores": collections.Counter(),
        "secondary_mismatches": collections.Counter(),
    }
        
    for line in view.stdout:
        # Parse the alignment JSON
        alignment = json.loads(line)
        
        # Count the alignment
        stats["total_reads"] += 1
        
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
                
    with open(stats_file, "w") as stats_handle:
        # Save the stats as JSON
        RealTimeLogger.get().info("Saving JSON stats: {}".format(stats))
        json.dump(stats, stats_handle)      
        
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
    
    RealTimeLogger.start_master(options)
    
    # Pre-read the input file so we don't try to send file handles over the
    # network.
    options.server_list = list(options.server_list)
    
    # Make a root job
    root_job = Job.wrapJobFn(run_all_alignments, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

