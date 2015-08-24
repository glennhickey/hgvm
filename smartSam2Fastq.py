#!/usr/bin/env python2.7
"""
smartSam2Fastq.py: turn sorted-by-names SAM input into a properly deduplicated 
FASTQ

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

    
def parse_and_deduplicate_sam(sam_input):
    """
    Given a source of input SAM lines, yield all the lines for unique reads.
    Parses lines into (name, flags, sequence, quality, mapping info) tuples,
    and un-reverses reads (and clears the reversed flag).
    
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
    
    # What was the name of the last read
    last_name = None
    
    # This holds a set of (forward sequence, is_read1, is_read2) tuples for a
    # given read name.
    found_reads = set()
    
    for line in sam_input:
        if line.startswith("@"):
            # Pass headers through
            yield line
            continue
    
        # Parse out the fields
        parts = line.split("\t")
        # Get the template name
        name = parts[0]
        
        if name != last_name:
            # Start a new state
            last_name = name
            found_reads = set()
        
        # Get the flags
        flags = int(parts[1])
        # Get the mapping info
        mapping_info = (parts[2], int(parts[3]), int(parts[4]))
        # Get the read sequence
        sequence = parts[9]
        # And the quality
        quality = parts[10]
        
        if flags & BAM_FREVERSE:
            # The sequence is given here as the reverse complement, so flip it
            # and the quality.
            sequence = reverse_complement(sequence)
            quality = quality[::-1]
            flags &= ~BAM_FREVERSE
            
        # Define the identity tuple
        identity = (sequence, flags & BAM_FREAD1, flags & BAM_FREAD2)
        if identity not in found_reads:
            # We don't have a copy of this yet
            found_reads.add(identity)
            
            yield (name, flags, sequence, quality, mapping_info)
            
def pair_up(records):
    """
    Pairs up deduplicated (name, flags, sequence, quality, mapping contig)
    records. Only yields matched pairs that come one after the other. Adds /1 to
    name for READ1 reads, and /2 for READ2 reads.
    
    Yields pairs of record lists
    
    """

    last_record = None
    
    for record in records:
    
        if last_record is not None and last_record[0] == record[0]:
            # We may have a matched pair
            
            if(record[1] & BAM_FREAD1 != last_record[1] & BAM_FREAD1 or
                record[1] & BAM_FREAD2 != last_record[1] & BAM_FREAD2):
                
                # One is a READ1 and one is a READ2
                
                # Make a mutable pair
                pair = [list(last_record), list(record)]
                
                for r in pair:
                    # Add the end tags to each name
                    if r[1] & BAM_FREAD1:
                        r[0] += "/1"
                    if r[1] & BAM_FREAD2:
                        r[0] += "/2"
                    
                # Yield them paired
                yield pair
                
                # Neither of these can be paired again
                last_record = None
                
            else:
                # Complain there's something wrong
                sys.stderr.write("Got improperly paired reads with name "
                    "{}\n".format(last_record[0]))
                sys.stderr.write("Is READ1: {} {}\n".format(
                    record[1] & BAM_FREAD1, last_record[1] & BAM_FREAD1))
                sys.stderr.write("Is READ2: {} {}\n".format(
                    record[1] & BAM_FREAD2, last_record[1] & BAM_FREAD2))
                sys.stderr.write("Read A:{}\n".format(record))
                sys.stderr.write("Read B:{}\n".format(last_record))
                
                # Just skip this read for now, until a proper partner comes
                # along for the other read, or some other name comes up.
                continue
                    
        else:
            # The last record doesn't match this one's name. It clearly has no
            # pair, so drop it. Remember this record instead.
            last_record = record
            
def write_fastq(stream, record):
    """
    Write the given record as FASTQ to the given stream
    """
    
    # Unpack and format the record
    (name, _, sequence, quality, _) = record
    stream.write("@{}\n{}\n+\n{}\n".format(name, sequence, quality))
    
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
    
    for (record1, record2) in pair_up(parse_and_deduplicate_sam(options.input_sam)):
        # Split up the records to their files
        write_fastq(options.fq1, record1)
        write_fastq(options.fq2, record2)
        
    
    
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

