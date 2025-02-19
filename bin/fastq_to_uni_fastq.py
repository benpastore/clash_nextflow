import sys
import gzip
import os

def fastq_to_uni(fastq, outfile) : 

    """
    Collapse PCR Dups in Fastq file output uniq fastq
    
    Input : fastq

    Output : uniq fastq
    """

    if ".gz" in fastq : 
        opener = gzip.open
        read_type = 'rt'
        print("reading in gz file")
    else : 
        opener = open
        read_type = 'r'
        print("reading uncompressed file")

    with opener(fastq, read_type) as f : 

        seqs = {}
        place = 0
        lines = []
        for line in f : 

            place += 1
            l = line.strip()
            lines.append(l)

            if place == 4 : 
                seq = lines[1]

                if seq in seqs.keys() : 
                    seqs[seq] += 1

                else : 
                    seqs[seq] = 1
                
                place = 0
                lines = []
    f.close()

    if os.path.exists(outfile) : 
        os.remove(outfile)

    f = open(outfile, 'a')

    for key, value in seqs.items() : 

        # sequence is key , count is value
        seq = key
        count = value
        qualscore = "I"*len(seq)

        f.write( "@{}:{}".format(seq, count) + "\n" + seq + "\n" + "+" + "\n" + qualscore + "\n")
    
    f.close() 

def main() : 

    fastq = sys.argv[1]
    outfile = sys.argv[2]

    fastq_to_uni(fastq, outfile)

if __name__ == "__main__" : 

    main()