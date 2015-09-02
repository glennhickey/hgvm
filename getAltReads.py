#!/usr/bin/env python2.7
"""
getAltReads.py: download reads for regions of interest from the 1000 Genomes FTP
server in parallel, using Toil.

"""

from toil.job import Job
import argparse, sys, os, os.path, random, collections, shutil, itertools, glob
import urllib2, urlparse, ftplib, fnmatch
import tsv

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
    parser.add_argument("--reference_metadata", 
        default="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCA_000001405.17_GRCh38.p2/"
        "GCA_000001405.17_GRCh38.p2_assembly_structure/"
        "all_alt_scaffold_placement.txt",
        help="URL to download the reference metadata from")
    parser.add_argument("--regions", nargs="*", 
        default=["BRCA1", "BRCA2", "CENX", "MHC", "SMA", "LRC_KIR"],
        help="region names to download reads for")
    parser.add_argument("--sample_ftp_root",
        default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data", 
        help="FTP directory to scan for samples")
    parser.add_argument("--sample_pattern", default="NA*", 
        help="fnmatch-style pattern for sample names")
    parser.add_argument("--sample_limit", type=int, default=100, 
        help="number of matching samples to download")
    parser.add_argument("out_dir",
        help="output directory to create and fill with per-region BAM files")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)



class HelloWorld(Job):
    def __init__(self):
        Job.__init__(self,  memory=100000, cores=2, disk=20000)
    def run(self, fileStore):
        fileId = fileStore.getEmptyFileStoreID()
        self.addChild(Job.wrapJobFn(childFn, fileId,
            cores=1, memory="1M", disk="10M"))
        self.addFollowOn(FollowOn(fileId))

def childFn(target, fileID):
    with target.fileStore.updateGlobalFileStream(fileID) as file:
        file.write("Hello, World!")

class FollowOn(Job):
    def __init__(self,fileId):
        Job.__init__(self)
        self.fileId=fileId
    def run(self, fileStore):
        tempDir = fileStore.getLocalTempDir()
        tempFilePath = "/".join([tempDir,"LocalCopy"])
        with fileStore.readGlobalFileStream(self.fileId) as globalFile:
            with open(tempFilePath, "w") as localFile:
               localFile.write(globalFile.read())

def downloadAllReads(job, options):
    """
    Download all the reads for the regions.
    
    """
    
    # First make the output directory
    if not os.path.exists(options.out_dir):
        try:
            # Make it if it doesn't exist
            os.makedirs(options.out_dir)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass
    
    # Whatever happens, it needs to exist here
    assert(os.path.exists(options.out_dir) and os.path.isdir(options.out_dir))
    
    # Holds the chromosome number for each region?
    region_chromosomes = {}
    # Holds the minimum start position for each region on its chromosome
    region_starts = collections.defaultdict(lambda: float("inf"))
    # Holds the maximum stop position for each region on its chromosome
    region_stops = collections.defaultdict(lambda: float("-inf"))
    
    # Holds the contig:start-end string for each alt in a given region
    # The reference range gets added in last
    ranges_by_region = collections.defaultdict(list)
    
    # Read the reference database
    database = tsv.TsvReader(urllib2.urlopen(options.reference_metadata))
    
    for parts in database:
        print parts
        # Parse out all the info for this alt and its parent chromosome
        region_name = parts[7]
        # Grab the chromosome ("1" or "X") that's the parent
        parent_chromosome = parts[5]
        parent_start = int(parts[11])
        parent_stop = int(parts[12])
        alt_contig = parts[3]
        alt_start = int(parts[9])
        alt_stop = int(parts[10])
        
        # Note the region start, stop, and parent chromosome number
        region_chromosomes[region_name] = parent_chromosome
        region_starts[region_name] = min(region_starts[region_name],
            parent_start)
        region_stops[region_name] = max(region_stops[region_name],
            parent_stop)
            
        # Turn the alt name into the proper format (GL000251.2 to
        # chr6_GL000251v2_alt)
        name_parts=alt_contig.split(".")
        fixed_alt_contig = "chr{}_{}v{}_alt".format(parent_chromosome,
            name_parts[0], name_parts[1])
            
        # Add it to the list for its region
        ranges_by_region[region_name].append("{}:{}-{}".format(
            fixed_alt_contig, alt_start, alt_stop))
                
    for region_name in region_chromosomes.iterkeys():
        # Add in the reference ranges that all the alts are alternatives to
        ranges_by_region[region_name].append("chr{}:{}-{}".format(
            region_chromosomes[region_name], region_starts[region_name], 
            region_stops[region_name]))
            

    # TODO: Download the sample list
    # Parse the sample directory to get the host and path
    ftp_info = urlparse.urlparse(options.sample_ftp_root)
    assert(ftp_info.scheme == "ftp")
    
    # Connect to the server
    ftp = ftplib.FTP(ftp_info.netloc, ftp_info.username, ftp_info.password)
    
    # Log in
    ftp.login()
        
    # Go to the right directory
    ftp.cwd(ftp_info.path)
    
    # Grab the right number of sample names that match the pattern from the
    # FTP directory listing.
    sample_names = itertools.islice(
        (n for n in ftp.nlst() if fnmatch.fnmatchcase(n,
        options.sample_pattern)), options.sample_limit)
    
    # Compose URLs for all of them and store them under the sample names
    sample_urls = {name: "{}/{}".format(options.sample_ftp_root, name)
        for name in sample_names}
    
    for region_name in region_chromosomes.iterkeys():
        for sample_name, sample_url in sample_urls.iteritems():
            
            # Where will this sample's BAM for this region go?
            bam_filename = "{}/{}/{}.bam".format(options.out_dir, region_name,
                sample_name)
            
            # Now kick off a job to download all the ranges for the region in
            # parallel for this sample, and then concatenate them together. Tell
            # it to save the results to a file on a shared filesystem.
            job.addChildFn(downloadRegion, region_name, sample_url,
                ranges_by_region[region_name], bam_filename, 
                cores=1, memory="4G", disk="50G")
   
def downloadRegion(job, region_name, sample_url, range_list, bam_filename):
    """
    Download all the ranges given for the given region from the given sample
    URL, concatenate them, and save them to the given BAM.
    
    """
    
    # TODO: implement
    print("Supposed to download {} ranges {} to {}".format(sample_url,
        range_list, bam_filename))

def main():
    options = parse_args(sys.argv) # This holds the nicely-parsed options object
    
    # Run Toil
    Job.Runner.startToil(Job.wrapJobFn(downloadAllReads, options,
        cores=1, memory="1G", disk=0),  options)

if __name__=="__main__":
    sys.exit(main())
