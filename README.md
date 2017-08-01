# fastmd
A Python script that computes missing "MD" tag values for SAM alignment files. Very fast when compared to "samtools calmd".

## The 'MD' tag in SAM files
The MD tag is an additional piece of information in a SAM alignment revealing the exact mismatches between the read and reference. Most programs that output SAM will include this MD tag automatically with the alignment for each read, while there are still some that unfortunately do not. The information in the MD tag is needed by a number of downstream analysis programs such as variant callers.

## usage
`fastmd.py reference.fasta input.sam output.sam`
