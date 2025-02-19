import pandas as pd
import subprocess
from multiprocessing import Pool
import re
from Bio import SeqIO  # For parsing FASTA files

# -------------------------------
# File paths
# -------------------------------
input_file = "/fs/ess/PCON0160/ben/pipelines/clash_nextflow/test.tsv"  # TSV file with gene coordinates
transcript_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/ce_ws279.linc.pseudo.pc.repbase.fa"  # Transcript sequences
smrna_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNA_reference.15ntextension.fa"  # smRNA sequences
output_file = "RNAduplex_results.tsv"

# -------------------------------
# (1) Load FASTA sequences into dictionaries
# -------------------------------
def load_fasta_sequences(fasta_file):
    """Parses a FASTA file and returns a dictionary {sequence_id: sequence}."""
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = str(record.seq).upper()  # Convert to uppercase
    return fasta_dict

transcript_sequences = load_fasta_sequences(transcript_fasta)
smrna_sequences = load_fasta_sequences(smrna_fasta)

# -------------------------------
# (2) Load the TSV file into a DataFrame
# -------------------------------
df = pd.read_csv(input_file, sep="\t")
# Uncomment the following line to restrict to the first 100 rows for testing:
df = df.head(100)

# -------------------------------
# Functions to process sequences
# -------------------------------
def extend_gene_mapped(row):
    """Extend the target sequence (gene_mapped) by ±5 nt using gene_start and gene_end."""
    gene_id = row["gene_id"]
    if gene_id in transcript_sequences:
        full_seq = transcript_sequences[gene_id]
        start = max(0, row["gene_start"] - 5)
        end = min(len(full_seq), row["gene_end"] + 5)
        return full_seq[start:end]
    return row["gene_mapped"]

def fix_smrna_sequence(row):
    """If smRNA_start is 1, 2, or 3, prepend the missing nucleotides from the full smRNA sequence."""
    smrna_id = row["smRNA_id"]
    smrna_start = row["smRNA_start"]
    if smrna_id in smrna_sequences:
        full_smrna_seq = smrna_sequences[smrna_id]
        if smrna_start in [1, 2, 3]:
            corrected_seq = full_smrna_seq[:smrna_start] + row["smRNA_seq"]
            return corrected_seq
    return row["smRNA_seq"]

# Process sequences: extend target and fix smRNA, then convert T->U.
df["extended_gene_mapped"] = df.apply(extend_gene_mapped, axis=1)
df["corrected_smRNA_seq"] = df.apply(fix_smrna_sequence, axis=1)
df["extended_gene_mapped"] = df["extended_gene_mapped"].str.replace("T", "U")
df["corrected_smRNA_seq"] = df["corrected_smRNA_seq"].str.replace("T", "U")

# Create the RNAduplex input in the format "target_seq&smrna_seq"
df["duplex_input"] = df["extended_gene_mapped"] + "&" + df["corrected_smRNA_seq"]

# Remove duplicate duplex inputs, if any.
df_unique = df.drop_duplicates(subset=["duplex_input"]).copy()

# -------------------------------
# Bulge Detection Function
# -------------------------------
def detect_bulges(dot_bracket):
    """
    Finds true bulges in a dot-bracket structure, distinguishing them from single-stranded regions.
    Returns a tuple: (target_bulges, smrna_bulges)
    """
    if "&" not in dot_bracket:
        return "Invalid Structure", "Invalid Structure"
    try:
        split_index = dot_bracket.index("&")
        target_structure = dot_bracket[:split_index]
        smrna_structure = dot_bracket[split_index+1:]
        min_length = min(len(target_structure), len(smrna_structure))
        def find_true_bulges(structure, complementary_structure, offset=0):
            bulge_positions = []
            for match in re.finditer(r"\.+", structure):
                start, end = match.start(), match.end() - 1
                if end >= min_length:
                    continue
                paired_bases = sum(1 for i in range(start, end+1) if i < min_length and complementary_structure[i] in "()")
                if paired_bases > 0:
                    bulge_positions.append(f"Bulge at {start+offset}-{end+offset}")
            return ", ".join(bulge_positions) if bulge_positions else "No Bulges"
        target_bulges = find_true_bulges(target_structure, smrna_structure)
        smrna_bulges = find_true_bulges(smrna_structure, target_structure, offset=0)
        return target_bulges, smrna_bulges
    except Exception as e:
        print(f"Error processing structure: {dot_bracket} - {e}")
        return "Error", "Error"

# -------------------------------
# RNAduplex Function
# -------------------------------
def run_rnaduplex(seq_pair):
    """
    Runs RNAduplex on a sequence pair.
    Expects input in the form "target_seq&smrna_seq".
    Returns a tuple:
      (seq_pair, Duplex Structure, Free Energy, Target Bulges, smrna Bulges)
    """
    try:
        # Split and strip sequences
        target_seq, smrna_seq = [s.strip() for s in seq_pair.split("&")]
        duplex_input = f"{target_seq}\n{smrna_seq}\n"
        
        result = subprocess.run(
            ["RNAduplex", "-T", "20"],
            input=duplex_input,
            capture_output=True,
            text=True
        )
        output = result.stdout.strip()
        if not output:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        output_lines = output.split("\n")
        if len(output_lines) == 1:
            line = output_lines[0]
            if " (" in line:
                structure, energy = line.split(" (")
                energy = energy.strip(") ")
            else:
                return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        elif len(output_lines) >= 2:
            # The second line is expected to have the dot-bracket and energy.
            structure_line = output_lines[1].strip()
            if " (" in structure_line:
                structure, energy = structure_line.split(" (")
                energy = energy.strip(") ")
            else:
                structure = structure_line
                energy = "NA"
        else:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        
        # Remove extra information from the structure.
        # Keep only the allowed characters: . ( ) &
        structure_clean = "".join(ch for ch in structure if ch in ".()&")
        
        target_bulges, smrna_bulges = detect_bulges(structure_clean)
        return seq_pair, structure_clean, energy, target_bulges, smrna_bulges
    except Exception as e:
        print(f"Error running RNAduplex for {seq_pair}: {e}")
        return seq_pair, "Exception", "Exception", "Exception", "Exception"

# -------------------------------
# Tail Base-Pairing Analysis Function
# -------------------------------
def is_tail_base_paired(row):
    """
    Checks if the non-templated 3' tail (from the 'tail' column)
    participates in base pairing in the RNAduplex structure.
    """
    tail_seq = row["tail"]  # For example, "U", "AA", or "*" (no tail)
    structure = row["Duplex Structure"]
    if not isinstance(structure, str) or "&" not in structure:
        return "Invalid Structure"
    if tail_seq == "*":
        return "No Tail"
    tail_length = len(tail_seq)
    split_index = structure.index("&")
    smrna_structure = structure[split_index+1:]
    if len(smrna_structure) < tail_length:
        return "Tail Out of Bounds"
    tail_structure = smrna_structure[-tail_length:]
    if "(" in tail_structure or ")" in tail_structure:
        return "Base-Paired"
    else:
        return "Unpaired"

# -------------------------------
# Run RNAduplex with multiprocessing
# -------------------------------
NUM_CORES = 8  # Adjust based on your system
with Pool(NUM_CORES) as p:
    results = p.map(run_rnaduplex, df_unique["duplex_input"].tolist())

# Create a DataFrame from RNAduplex results
df_results = pd.DataFrame(results, columns=["duplex_input", "Duplex Structure", "Free Energy", "Target Bulges", "smRNA Bulges"])

# Merge RNAduplex results with the original DataFrame
df_final = df.merge(df_results, on="duplex_input", how="left")

# Add Tail Base-Paired Analysis
df_final["Tail Base-Paired?"] = df_final.apply(is_tail_base_paired, axis=1)
print(df_final[["tail", "Duplex Structure", "Tail Base-Paired?"]].head())

# Save final results to TSV
df_final.to_csv(output_file, sep="\t", index=False)
print(f"RNAduplex predictions complete! Results saved to {output_file}")
