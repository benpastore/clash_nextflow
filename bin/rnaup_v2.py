import pandas as pd
import subprocess
from multiprocessing import Pool
import re
from Bio import SeqIO  # For parsing FASTA files

# Load input files
input_file = "/fs/ess/PCON0160/ben/pipelines/clash_nextflow/test.tsv"  # TSV file with gene coordinates
transcript_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/ce_ws279.linc.pseudo.pc.repbase.fa"  # FASTA file with transcript sequences
smrna_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNA_reference.15ntextension.fa"  # FASTA file with smRNA sequences
output_file = "RNAcofold_results_with_bulges.tsv"

# ✅ (1) Function to Load Sequences from FASTA into a Dictionary
def load_fasta_sequences(fasta_file):
    """Parses a FASTA file and returns a dictionary {sequence_id: sequence}."""
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = str(record.seq).upper()  # Convert to uppercase
    return fasta_dict

# Load transcript sequences and smRNA sequences
transcript_sequences = load_fasta_sequences(transcript_fasta)
smrna_sequences = load_fasta_sequences(smrna_fasta)

# ✅ (2) Load TSV file into a pandas DataFrame
df = pd.read_csv(input_file, sep="\t")
df = df.head(100)

# Function to extend the target site (gene_mapped) using gene_start and gene_end
def extend_gene_mapped(row):
    gene_id = row["gene_id"]
    if gene_id in transcript_sequences:
        full_seq = transcript_sequences[gene_id]
        start = max(0, row["gene_start"] - 5)  # Ensure start doesn't go below 0
        end = min(len(full_seq), row["gene_end"] + 5)  # Ensure end doesn't exceed transcript length
        return full_seq[start:end]
    return row["gene_mapped"]  # Default to original if transcript not found

# Function to fix smRNA sequences using smRNA_start
def fix_smrna_sequence(row):
    smrna_id = row["smRNA_id"]
    smrna_start = row["smRNA_start"]
    
    if smrna_id in smrna_sequences:
        full_smrna_seq = smrna_sequences[smrna_id]  # Get full sequence from FASTA
        if smrna_start in [1, 2, 3]:  # If smRNA starts at 1, 2, or 3
            corrected_seq = full_smrna_seq[:smrna_start] + row["smRNA_seq"]  # Add missing bases
            return corrected_seq
    return row["smRNA_seq"]  # Default to original if no correction needed

# Extend gene_mapped sequences by ±5 nt
df["extended_gene_mapped"] = df.apply(extend_gene_mapped, axis=1)

# Correct smRNA sequences based on smRNA_start
df["corrected_smRNA_seq"] = df.apply(fix_smrna_sequence, axis=1)

# Convert T to U in both extended_gene_mapped and corrected_smRNA_seq
df["extended_gene_mapped"] = df["extended_gene_mapped"].str.replace("T", "U")
df["corrected_smRNA_seq"] = df["corrected_smRNA_seq"].str.replace("T", "U")

# Create RNAcofold input format: "extended_gene_mapped&corrected_smRNA_seq"
df["cofold_input"] = df["extended_gene_mapped"] + "&" + df["corrected_smRNA_seq"]

# ✅ Remove duplicates based on extended_gene_mapped and corrected_smRNA_seq
df_unique = df.drop_duplicates(subset=["cofold_input"]).copy()

def detect_bulges(dot_bracket):
    """Finds bulges (multiple unpaired bases) in dot-bracket notation separately for the target and smRNA."""
    if "&" not in dot_bracket:
        return "Invalid Structure", "Invalid Structure"

    # Identify the split position
    split_index = dot_bracket.index("&")  # The "&" separates target RNA and smRNA
    target_structure = dot_bracket[:split_index]  # Dot-bracket notation of target RNA
    smrna_structure = dot_bracket[split_index + 1:]  # Dot-bracket notation of smRNA

    # Function to find bulges in a given structure
    def find_bulges(structure, offset=0):
        bulge_positions = []
        for match in re.finditer(r"\.+", structure):  # Find consecutive dots (unpaired bases)
            if len(match.group()) > 1:  # Bulges must be >1 nt long
                start = match.start() + offset  # Adjust index
                end = match.end() - 1 + offset
                bulge_positions.append(f"Bulge at {start}-{end}")
        return ", ".join(bulge_positions) if bulge_positions else "No Bulges"

    # Detect bulges separately
    target_bulges = find_bulges(target_structure)
    smrna_bulges = find_bulges(smrna_structure, offset=0)  # smRNA starts at 0

    return target_bulges, smrna_bulges

def run_rnacofold(seq_pair):
    """Runs RNAcofold on a sequence pair and extracts the folding pattern and free energy."""
    try:
        # Debugging step: Print input sequence before running RNAcofold
        #print(f"Running RNAcofold for: {seq_pair}")

        # Run RNAcofold
        result = subprocess.run(
            ["RNAcofold", "-T", "20"],  # Run at 20°C
            input=seq_pair,
            capture_output=True,
            text=True
        )

        # Debugging step: Print output from RNAcofold
        #print(f"RNAcofold Output:\n{result.stdout}\nErrors:\n{result.stderr}")

        # Parse RNAcofold output
        output_lines = result.stdout.strip().split("\n")
        if len(output_lines) < 2 or "ERROR" in result.stderr:
            return seq_pair, "RNAcofold Failed", "RNAcofold Failed", "No Structure", "No Bulges"

        # Extract folding structure and free energy
        structure_line = output_lines[1]
        structure, energy = structure_line.split(" (")
        energy = energy.strip(") kcal/mol")

        # Detect bulges separately for the target and smRNA
        target_bulges, smrna_bulges = detect_bulges(structure.strip())

        return seq_pair, structure.strip(), energy.strip(), target_bulges, smrna_bulges
    
    except Exception as e:
        print(f"Error running RNAcofold for {seq_pair}: {e}")
        return seq_pair, "Exception", "Exception", "Exception", "Exception"


def generate_base_pairing_alignment(row):
    """Creates a visual representation of base pairing between target RNA and smRNA."""
    structure = row["Folding Pattern"]
    target_seq = row["extended_gene_mapped"]
    smrna_seq = row["corrected_smRNA_seq"]

    # Check for errors
    if "Error" in structure or "&" not in structure:
        return "Invalid Structure"

    # Extract target and smRNA structures
    split_index = structure.index("&")
    target_structure = structure[:split_index]
    smrna_structure = structure[split_index + 1:]

    # Align sequences
    paired_target = []
    alignment = []
    paired_smrna = []

    # Reverse smRNA to match the pairing orientation
    smrna_seq = smrna_seq[::-1]  # Reverse small RNA for 3' -> 5' alignment

    # Iterate over structure and create alignment
    smrna_index = 0
    for i, (t_base, t_struct) in enumerate(zip(target_seq, target_structure)):
        if t_struct in "().":  # Only process actual bases
            paired_target.append(t_base)
            if t_struct == "(":  # Paired base
                if smrna_index < len(smrna_seq):
                    s_base = smrna_seq[smrna_index]
                    if (t_base, s_base) in [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")]:
                        alignment.append("|")  # Watson-Crick pair
                    elif (t_base, s_base) in [("G", "U"), ("U", "G")]:
                        alignment.append("+")  # G-U wobble
                    else:
                        alignment.append(".")  # Mismatch
                    paired_smrna.append(s_base)
                    smrna_index += 1
                else:
                    alignment.append(".")
                    paired_smrna.append("-")  # No matching smRNA base
            else:
                alignment.append(".")  # Unpaired
                paired_smrna.append("-")  # No matching smRNA base

    # Reverse small RNA back for correct display
    paired_smrna = paired_smrna[::-1]

    # Create formatted alignment
    formatted_alignment = f"Target:    {''.join(paired_target)}\n" \
                          f"           {''.join(alignment)}\n" \
                          f"smRNA:     {''.join(paired_smrna)}"

    return formatted_alignment


# ✅ Function to Check if Tail is Base-Paired
def is_tail_base_paired(row):
    """Checks if the non-templated 3' nt (tail) participates in base pairing."""
    tail_seq = row["tail"]  # The tail sequence (e.g., "U", "AA", or "*")
    structure = row["Folding Pattern"]

    # Check if structure is valid and contains "&"
    if not isinstance(structure, str) or "&" not in structure:
        return "Invalid Structure"

    if tail_seq == "*":  # No tail present
        return "No Tail"

    # Get the length of the tail
    tail_length = len(tail_seq)
    smrna_length = len(row["corrected_smRNA_seq"])  # Full corrected smRNA sequence

    # Find the split index in dot-bracket notation
    split_index = structure.index("&")
    smrna_structure = structure[split_index + 1:]  # smRNA dot-bracket structure

    # Ensure we have enough bases to analyze
    if len(smrna_structure) < tail_length:
        return "Tail Out of Bounds"

    # Identify the dot-bracket symbols for the last 'tail_length' nucleotides
    tail_structure = smrna_structure[-tail_length:]

    # Check if tail is paired or unpaired
    if "(" in tail_structure or ")" in tail_structure:
        return "Base-Paired"
    else:
        return "Unpaired"

# Use multiprocessing to speed up RNAcofold runs
NUM_CORES = 8  # Adjust based on your system
with Pool(NUM_CORES) as p:
    results = p.map(run_rnacofold, df_unique["cofold_input"].tolist())

# Convert results to DataFrame
df_results = pd.DataFrame(results, columns=["cofold_input", "Folding Pattern", "Free Energy", "Target Bulges", "smRNA Bulges"])

# ✅ Merge back RNAcofold results with the original dataframe
df_final = df.merge(df_results, on="cofold_input", how="left")

# ✅ Add Tail Base-Pairing Analysis
df_final["Tail Base-Paired?"] = df_final.apply(is_tail_base_paired, axis=1)
print(df_final[["tail", "Folding Pattern", "Tail Base-Paired?"]].head())  # Debugging step

df_final["Base-Pairing Alignment"] = df_final.apply(generate_base_pairing_alignment, axis=1)

# Save to a new TSV file
df_final.to_csv(output_file, sep="\t", index=False)

print(f"RNAcofold predictions complete! Results saved to {output_file}")
