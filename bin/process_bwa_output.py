import pysam
import sys

# Open the SAM file
samfile = pysam.AlignmentFile(sys.argv[1], "r")
outfile = sys.argv[2]

# Open the output BED file for writing.
# In this example, we include extra columns for the aligned (target) and unaligned (piRNA) sequences.
with open(f"{outfile}", "w") as bedfile:
    for read in samfile:
        # Skip unmapped reads.
        if read.is_unmapped:
            continue

        # Get the read sequence
        read_seq = read.query_sequence

        # Initialize empty strings for the aligned (target) and unaligned (piRNA) portions.
        aligned_seq = ""
        unaligned_seq = ""
        pos = 0  # pointer into the read sequence

        # Loop through the CIGAR operations once.
        # CIGAR op codes of interest: 0 = match/mismatch (aligned), 4 = soft clip (unaligned)
        for op, length in read.cigartuples:
            if op == 4:
                # Soft-clipped segment
                unaligned_seq += read_seq[pos:pos+length]
                pos += length
            elif op == 0:
                # Matched (aligned) segment
                aligned_seq += read_seq[pos:pos+length]
                pos += length
            else:
                # For any other CIGAR operations, just advance the pointer.
                pos += length

        # For reporting, we assume that the aligned (M) portion corresponds to the target transcript.
        # (If you want to flag orientation differences, you could inspect the order of the operations.)

        # Extract mapping information for BED fields.
        # BED fields: chrom, chromStart, chromEnd, name, score, strand.
        chrom = read.reference_name
        chromStart = read.reference_start  # 0-based start
        chromEnd = read.reference_end      # end coordinate (non-inclusive)
        name = read.query_name
        score = read.mapping_quality
        strand = "-" if read.is_reverse else "+"

        # Write a BED line.
        # We include two extra columns: aligned_seq (target transcript) and unaligned_seq (piRNA).
        if strand == "+" : 
            bed_line = (
                f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{score}\t{strand}\t"
                f"{aligned_seq}\t{unaligned_seq}\n"
            )
            bedfile.write(bed_line)

samfile.close()
