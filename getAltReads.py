#!/usr/bin/env python2.7
"""
getAltReads.py: download reads for regions of interest from the 1000 Genomes FTP
server in parallel, using Toil.

"""

from toil.job import Job
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob

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
    parser.add_argument("--server_list", type=argparse.FileType("r"),
        help="TSV file continaing <region>\t<url> lines for servers to test")
    
    
    
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

def main():
    options = parse_args(sys.argv) # This holds the nicely-parsed options object
    
    # Run Toil
    Job.Runner.startToil(HelloWorld(),  options)

if __name__=="__main__":
    sys.exit(main())
