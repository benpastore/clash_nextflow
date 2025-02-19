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
# Optionally, restrict to a subset for testing:
# df = df.head(100)

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

# Process sequences and convert T -> U.
df["extended_gene_mapped"] = df.apply(extend_gene_mapped, axis=1)
df["corrected_smRNA_seq"] = df.apply(fix_smrna_sequence, axis=1)
df["extended_gene_mapped"] = df["extended_gene_mapped"].str.replace("T", "U")
df["corrected_smRNA_seq"] = df["corrected_smRNA_seq"].str.replace("T", "U")

# Create RNAduplex input in the form: "extended_gene_mapped&corrected_smRNA_seq"
df["duplex_input"] = df["extended_gene_mapped"] + "&" + df["corrected_smRNA_seq"]

# Remove duplicate duplex inputs.
df_unique = df.drop_duplicates(subset=["duplex_input"]).copy()

# -------------------------------
# Bulge detection function (updated)
# -------------------------------
def detect_bulges(dot_bracket):
    """
    Detects "true bulges" in a dot‐bracket structure and distinguishes them from simply unpaired terminal regions.
    
    A region in a given strand is reported as a bulge if:
      (1) It is internal in that strand (i.e. it does not extend to the very start or end).
      (2) It is flanked on both sides by paired bases (i.e. the base immediately before and immediately after the unpaired region are "(" or ")").
      (3) When comparing the unpaired region (or its overlapping portion, if the region is longer than the complementary strand) 
          with the complementary strand, at least one base is paired in that overlapping region.
          
    The function splits the structure at the "&" into target and smRNA parts and returns a tuple:
         (target_bulges, smrna_bulges).
    If no true bulges are found in a strand, it returns "No Bulges" for that strand.
    """
    if "&" not in dot_bracket:
        return "Invalid Structure", "Invalid Structure"
    try:
        # Split into target and smRNA parts.
        split_index = dot_bracket.index("&")
        target_structure = dot_bracket[:split_index]
        smrna_structure = dot_bracket[split_index + 1:]
        
        # Helper function: find true bulges in one strand.
        def find_true_bulges(structure, complementary_structure, offset=0):
            bulge_positions = []
            L = len(structure)
            L_comp = len(complementary_structure)
            for match in re.finditer(r"\.+", structure):
                start, end = match.start(), match.end() - 1  # end is inclusive
                # (1) Ignore stretches that touch the boundaries.
                if start == 0 or end == L - 1:
                    continue
                # (2) Require that the stretch is flanked by paired bases.
                if structure[start - 1] not in "()" or structure[end + 1] not in "()":
                    continue
                # (3) Determine the overlapping region in the complementary strand.
                # If the unpaired region extends beyond the complementary strand, only consider the overlapping portion.
                effective_start = start
                effective_end = min(end, L_comp - 1)
                # If there is no overlap, skip.
                if effective_end < effective_start:
                    continue
                overlap_length = effective_end - effective_start + 1
                # Count the number of paired bases in the complementary strand over the overlapping region.
                paired_bases = sum(1 for i in range(effective_start, effective_end + 1)
                                   if complementary_structure[i] in "()")
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
# RNAduplex function (with full-length padding for both strands)
# -------------------------------
def run_rnaduplex(seq_pair):
    """
    Runs RNAduplex on a sequence pair.
    Expects input in the form "target_seq&smrna_seq".
    Returns a tuple:
      (seq_pair, Final Duplex Structure, Free Energy, Target Bulges, smrna Bulges)
    
    The final duplex structure is modified so that:
      - The target portion is padded with dots to match the full length of the extended target sequence.
      - The smRNA portion is padded with dots to match the full length of the corrected smRNA sequence.
    """
    try:
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
            header_line = output_lines[0].strip()  # Unused in this version.
            structure_line = output_lines[1].strip()
            if " (" in structure_line:
                structure, energy = structure_line.split(" (")
                energy = energy.strip(") ")
            else:
                structure = structure_line
                energy = "NA"
        else:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        
        # Keep only allowed characters.
        structure_clean = "".join(ch for ch in structure if ch in ".()&")
        if "&" not in structure_clean:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        predicted_target_structure, predicted_smrna_structure = structure_clean.split("&")
        
        # Pad the target structure to match the full length of target_seq.
        full_target_len = len(target_seq)
        predicted_target_len = len(predicted_target_structure)
        pad_total_target = full_target_len - predicted_target_len
        if pad_total_target < 0:
            pad_total_target = 0
        pad_left_target = pad_total_target // 2
        pad_right_target = pad_total_target - pad_left_target
        padded_target_structure = ("." * pad_left_target) + predicted_target_structure + ("." * pad_right_target)
        
        # Pad the smRNA structure to match the full length of smrna_seq.
        full_smrna_len = len(smrna_seq)
        predicted_smrna_len = len(predicted_smrna_structure)
        pad_total_smrna = full_smrna_len - predicted_smrna_len
        if pad_total_smrna < 0:
            pad_total_smrna = 0
        pad_left_smrna = pad_total_smrna // 2
        pad_right_smrna = pad_total_smrna - pad_left_smrna
        padded_smrna_structure = ("." * pad_left_smrna) + predicted_smrna_structure + ("." * pad_right_smrna)
        
        final_structure = padded_target_structure + "&" + padded_smrna_structure
        final_structure_clean = "".join(ch for ch in final_structure if ch in ".()&")
        target_bulges, smrna_bulges = detect_bulges(final_structure_clean)
        return seq_pair, final_structure_clean, energy, target_bulges, smrna_bulges
    except Exception as e:
        print(f"Error running RNAduplex for {seq_pair}: {e}")
        return seq_pair, "Exception", "Exception", "Exception", "Exception"

# -------------------------------
# Tail base-pairing analysis function
# -------------------------------
def is_tail_base_paired(row):
    """
    Checks if the non-templated 3' tail (from the 'tail' column)
    participates in base pairing in the RNAduplex structure.
    """
    tail_seq = row["tail"]
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
# Function to extract target interaction region
# -------------------------------
def extract_target_interaction_region(row, extension=2):
    """
    Extracts the region of the target sequence (extended_gene_mapped) that is predicted to interact with the smRNA.
    
    The function examines the target portion of the Duplex Structure (i.e. the part before the "&"),
    finds the first and last indices where the structure is not ".", and then extends those boundaries by
    'extension' nucleotides (if possible). It then returns the corresponding substring from the target sequence.
    
    This region should represent the predicted interacting region on the target that can be annotated.
    """
    target_seq = row["extended_gene_mapped"]
    duplex_structure = row["Duplex Structure"]
    if not isinstance(duplex_structure, str) or "&" not in duplex_structure:
        return ""
    target_structure = duplex_structure.split("&")[0]
    indices = [i for i, ch in enumerate(target_structure) if ch != '.']
    if not indices:
        return ""
    start_idx = min(indices)
    end_idx = max(indices)
    # Extend boundaries by 'extension' nucleotides if available.
    start_idx = max(0, start_idx - extension)
    end_idx = min(len(target_seq) - 1, end_idx + extension)
    return target_seq[start_idx:end_idx+1]

# -------------------------------
# Run RNAduplex using multiprocessing
# -------------------------------
NUM_CORES = 8  # Adjust as needed
with Pool(NUM_CORES) as p:
    results = p.map(run_rnaduplex, df_unique["duplex_input"].tolist())

# Create a DataFrame from RNAduplex results.
df_results = pd.DataFrame(results, columns=["duplex_input", "Duplex Structure", "Free Energy", "Target Bulges", "smRNA Bulges"])

# Merge RNAduplex results with the original DataFrame.
df_final = df.merge(df_results, on="duplex_input", how="left")

# Add Tail Base-Paired Analysis.
df_final["Tail Base-Paired?"] = df_final.apply(is_tail_base_paired, axis=1)

# Add the Target Interaction Region (for genome annotation).
df_final["Target Interaction Region"] = df_final.apply(extract_target_interaction_region, axis=1)

# Debug print (optional)
print(df_final[["extended_gene_mapped", "Duplex Structure", "Target Interaction Region"]].head())

# Save the final results to a TSV file.
df_final.to_csv(output_file, sep="\t", index=False)
print(f"RNAduplex predictions complete! Results saved to {output_file}")
