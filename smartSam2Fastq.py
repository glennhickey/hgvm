#!/usr/bin/env python2.7
"""
smartSam2Fastq.py: turn sorted-by-names SAM input into a properly deduplicated
FASTQ. Also work around the bug in bwa mem where some even-length alignments to
reverse strands of alts will have incorrect bases for one of the middle two
bases.

accounts for both secondary and supplementary alignments of the same read

"""

import argparse, sys, os, os.path, random, itertools, string
import doctest

# We need to do reverse complements
RC_TABLE=string.maketrans("ACGT", "TGCA")

def reverse_complement(dna):
    # Simplest way to do it according to
    # <http://stackoverflow.com/a/26615937/402891>
    return dna.translate(RC_TABLE)[::-1]

# Global SAM constants
BAM_FREVERSE=16
BAM_FREAD1=64
BAM_FREAD2=128

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
    
    # General options
    parser.add_argument("--input_sam", type=argparse.FileType("r"),
        default=sys.stdin,
        help="input SAM in name-sorted order")
    parser.add_argument("--fq1", type=argparse.FileType("w"),
        default=sys.stdout,
        help="FASTQ file to save all the READ1 reads in")
    parser.add_argument("--fq2", type=argparse.FileType("w"),
        default=sys.stdout,
        help="FASTQ file to save all the READ2 reads in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

class Read(object):
    """
    Represent a Read as reconstructed from an alignment.
    
    """
    def __init__(self, sam_line):
        """
        Parse the given SAM line and construct a read.
        
        """
        
        # Save the line
        self.line = sam_line
        
        # Parse out the fields
        parts = sam_line.split("\t")
        # Get the template name
        self.template = parts[0]
        
        # Grab the flags
        self.flags = int(parts[1])
        
        # What end are we (1, 2, or 0 for unpaired)
        if self.flags & BAM_FREAD1:
            if self.flags & BAM_FREAD2:
                # This shouldn't happen
                raise RuntimeError("Alignment flagged as both READ1 and READ2")
            
            # Otherwise we're READ1
            self.end = 1
        elif self.flags & BAM_FREAD2:
            # We're READ2
            self.end = 2
        else:
            # We're unpaired
            self.end = 0
            
        # Grab sequence and qualities
        self.sequence = parts[9]
        self.qualities = parts[10]
        
        
        
        if self.flags & BAM_FREVERSE:
            # Flip to the other strand by RCing sequence and reversing qualities
            self.sequence = reverse_complement(self.sequence)
            self.qualities = self.qualities[::-1]
            self.is_reverse = True
        else:
            self.is_reverse = False
            
        # Grab the contig we mapped to
        self.contig = parts[2]
        
        # Figure out if we are suspected to be corrupted by the bwa mem bug or
        # not. It seems to affect all alignments to alts, not just reverse strand
        # ones.
        self.is_suspect = ((True or self.flags & BAM_FREVERSE) and 
            self.contig.endswith("_alt") and
            len(self.sequence) % 2 == 0)
            
    def get_name(self):
        """
        Produce a name for the read based on the template and the end.
        
        """
        
        if self.end == 0:
            # Unpaired reads are named for the template
            return self.template
        else:
            # Paired reads get /1 and /2
            return "{}/{}".format(self.template, self.end)
            
    def __eq__(self, other):
        """
        Two Reads are equal when they are the same read, but other meta-info
        about the original alignment may differ.
        
        """
        
        if not isinstance(other, self.__class__):
            return False
        
        
        if self.template != other.template:
            return False
        if self.sequence != other.sequence:
            return False
        if self.qualities != other.qualities:
            return Fasle
        if self.end != other.end:
            return False
            
        return True
        
    def __ne__(self, other):
        return not (self.__eq__(other))
        
    def __gt__(self, other):
        """
        A Read is greater than another read if it should replace the other read
        (they are the same end of the same template, and the other read is
        suspect while this one is not, or this one has a longer sequence and
        isn't suspect).
        
        """
        
        return (self.template == other.template and self.end == other.end and
            (not self.is_suspect) and 
            (len(self.sequence) > len(other.sequence) or other.is_suspect))
            
    def __str__(self):
        """
        Turn this Read into a string for displaying.
        
        """
        
        mark = "?" if self.is_suspect else ""
        return "{} end {} on {}: {}{}".format(self.template, self.end,
            self.contig, self.sequence, mark)
        
def parse_and_deduplicate_sam(sam_input):
    """
    Given a source of input SAM lines, parses lines into Read objects, and
    deduplicates them, discarding suspect ones when non-suspect ones are
    available.
    
    Yields dicts form end number to Read object for each template.
    
    TODO: fix this test to not be ugly.
    >>> lines_in = ["ERR894727.320\\t2147\\tchr6_GL000254v2_alt\\t198561\\t60"
    ... "\\t126M"
    ... "\\t=\\t198968\\t533\\tNCACCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCTGTCTCAA"
    ... "AAAA"
    ... "CAAACAAACAAACAACAACAACAACAGAAAACAGGGTGCAGCCCACTCCTCCAGCACCTTGAATCTGGTG"
    ... "GGCT\\t'7<<BB<BBBB<B<BBB0<BBBBBB<<BBBBBBB7B<<BBB<0<BB<BB0<B<BBB<BBB<B<"
    ... "<"
    ... "<BB<BB<BB<BB<BBBBBBBBB<<0070BBBBB<B00<B<B0B07<00<B00B00<0<<<<<",
    ... "ERR894727.320\\t2147\\tchr6_GL000255v2_alt\\t198391\\t60\\t126M\\t=\\t"
    ... "198798"
    ... "\\t533\\tNCACCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCTGTCTCAAAAAACAAACA"
    ... "AACAAACAACAACAACAACAGAAAACAGGGTGCAGCCCACTCCTCCAGCACCTTGAATCTGGTGGGCT"
    ... "\\t"
    ... "'7<<BB<BBBB<B<BBB0<BBBBBB<<BBBBBBB7B<<BBB<0<BB<BB0<B<BBB<BBB<B<<<BB<BB"
    ... "<BB<BB<BBBBBBBBB<<0070BBBBB<B00<B<B0B07<00<B00B00<0<<<<<"]
    >>> len(list(parse_and_deduplicate_sam(lines_in)))
    1
    
    """
    
    # What was the template for the last read
    last_template = None
    
    # For this template, we keep the best read for each end we find.
    reads_by_end = {}
    
    for line in sam_input:
        if line.startswith("@"):
            # Pass headers through
            yield line
            continue
    
        # Parse the read
        read = Read(line)
        
        if read.template != last_template:
            
            all_ok = True
            for end_read in reads_by_end.itervalues():
                # Yield out the last state's reads.
                
                if(end_read.is_suspect):
                    sys.stderr.write("Only suspect alignments found for end "
                        "{} of template {}. Skipping.\n".format(end_read.end,
                        end_read.template))
                all_ok = False
              
            if all_ok:  
                yield reads_by_end
        
            # Start a new state
            last_template = read.template
            reads_by_end = {}
            
        
        if not reads_by_end.has_key(read.end):
            # This is the only read for this end so far
            reads_by_end[read.end] = read
        else:
            if reads_by_end[read.end] < read:
                # Replace the existing read
                reads_by_end[read.end] = read
            elif (not read.is_suspect and 
                len(read.sequence) >= len(reads_by_end[read.end].sequence) and 
                read != reads_by_end[read.end]):
                # We aren't suspect, we differ, and we can't replace the other
                # read.
                raise RuntimeError("Non-suspect alignments don't agree on end "
                    "{} of template {}:\n{}\n{}".format(read.end, read.template,
                    read.line, reads_by_end[read.end].line))
    
    all_ok = True
    for end_read in reads_by_end.itervalues():
        # Yield out the last state's reads.
        
        if(end_read.is_suspect):
            sys.stderr.write("Only suspect alignments found for end "
                "{} of template {}. Skipping.\n".format(end_read.end,
                end_read.template))
        all_ok = False
      
    if all_ok:  
        yield reads_by_end
            
def write_fastq(stream, read):
    """
    Write the given record as FASTQ to the given stream
    """
    
    # Unpack and format the record
    stream.write("@{}\n{}\n+\n{}\n".format(read.get_name(), read.sequence,
        read.qualities))
    
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
    
    for reads_by_end in parse_and_deduplicate_sam(options.input_sam):
        if not (reads_by_end.has_key(1) and reads_by_end.has_key(2)):
            # Skip unpaired reads
            continue
            
        # Split up the reads to their files
        write_fastq(options.fq1, reads_by_end[1])
        write_fastq(options.fq2, reads_by_end[2])
        
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

