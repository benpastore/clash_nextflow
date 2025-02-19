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
    """
    Fixes the smRNA sequence by prepending the missing nucleotides if smRNA_start is 1, 2, or 3.
    After fixing, it checks whether the first nucleotide is T (or U).
    If not, it examines the first three nucleotides of the full smRNA sequence
    and replaces the first nucleotide with the first one that is T or U.
    """
    smrna_id = row["smRNA_id"]
    smrna_start = row["smRNA_start"]
    if smrna_id in smrna_sequences:
        full_smrna_seq = smrna_sequences[smrna_id]
        # If smRNA_start is 1, 2, or 3, prepend the missing bases.
        if smrna_start > 0 :
            corrected_seq = full_smrna_seq[:smrna_start] + row["smRNA_seq"]
        else:
            corrected_seq = row["smRNA_seq"]
        # Check if the corrected sequence starts with T or U.
        # (Note: After converting, T will become U, but we check before conversion.)
        if corrected_seq[0] not in ["T", "U"]:
            # Look at the first three nucleotides of the full smRNA sequence.
            for i in range(3):
                if full_smrna_seq[i] in ["T", "U"]:
                    corrected_seq = full_smrna_seq[i] + corrected_seq[1:]
                    break
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
# Bulge detection function (unchanged from previous update)
# -------------------------------
def detect_bulges(dot_bracket):
    """
    Detects "true bulges" in a dot bracket structure and distinguishes them from terminal unpaired regions.
    
    A region is reported as a bulge if:
      (1) It is internal (i.e. does not extend to the very start or end).
      (2) It is flanked on both sides by paired bases.
      (3) In the overlapping region of the complementary strand, at least one base is paired.
      (4) The unpaired stretch is at least 2 nucleotides long.
      
    Returns a tuple: (target_bulges, smrna_bulges)
    """

    if "&" not in dot_bracket:
        return "Invalid Structure", "Invalid Structure"
    try:
        split_index = dot_bracket.index("&")
        target_structure = dot_bracket[:split_index]
        smrna_structure = dot_bracket[split_index+1:]
        def find_true_bulges(structure, complementary_structure, offset=0):
            bulge_positions = []
            L = len(structure)
            L_comp = len(complementary_structure)
            for match in re.finditer(r"\.+", structure):
                region_length = match.end() - match.start()
                if region_length < 2:
                    continue
                start, end = match.start(), match.end() - 1  # end inclusive
                if start == 0 or end == L - 1:
                    continue
                if structure[start - 1] not in "()" or structure[end + 1] not in "()":
                    continue
                effective_start = start
                effective_end = min(end, L_comp - 1)
                if effective_end < effective_start:
                    continue
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
# RNAduplex function with full-length padding for both strands
# -------------------------------
def run_rnaduplex(seq_pair):
    """
    Runs RNAduplex on a sequence pair.
    Expects input in the form "target_seq&smrna_seq".
    Returns a tuple:
      (seq_pair, Final Duplex Structure, Free Energy, Target Bulges, smrna Bulges)
    
    The final duplex structure is modified so that:
      - The target portion is padded with dots to match the full length of target_seq.
      - The smRNA portion is padded with dots to match the full length of smrna_seq.
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
            header_line = output_lines[0].strip()  # Not used here.
            structure_line = output_lines[1].strip()
            if " (" in structure_line:
                structure, energy = structure_line.split(" (")
                energy = energy.strip(") ")
            else:
                structure = structure_line
                energy = "NA"
        else:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        
        # Clean structure: keep only allowed characters.
        structure_clean = "".join(ch for ch in structure if ch in ".()&")
        if "&" not in structure_clean:
            return seq_pair, "RNAduplex Failed", "RNAduplex Failed", "No Bulges", "No Bulges"
        predicted_target_structure, predicted_smrna_structure = structure_clean.split("&")
        
        # Pad the target structure.
        full_target_len = len(target_seq)
        predicted_target_len = len(predicted_target_structure)
        pad_total_target = full_target_len - predicted_target_len
        if pad_total_target < 0:
            pad_total_target = 0
        pad_left_target = pad_total_target // 2
        pad_right_target = pad_total_target - pad_left_target
        padded_target_structure = ("." * pad_left_target) + predicted_target_structure + ("." * pad_right_target)
        
        # Pad the smRNA structure.
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
# New Function: Annotate Base-Pair Counts in smRNA Regions
# -------------------------------
def annotate_bp_regions(duplex_structure):
    """
    Given a duplex structure in the form "target_structure&smrna_structure",
    this function extracts the smRNA part and counts the number of base-paired
    nucleotides (i.e. "(" or ")") in three regions:
      - Seed region: nucleotides 2-8 (i.e. indices 1 to 7)
      - Central region: nucleotides 9-13 (i.e. indices 8 to 12)
      - Supplemental region: nucleotides 14-end (i.e. index 13 onward)
    Returns a tuple: (seed_bp, central_bp, supplemental_bp)
    """
    if not isinstance(duplex_structure, str) or "&" not in duplex_structure:
        return None, None, None
    try:
        # Split the structure into target and smRNA parts.
        parts = duplex_structure.split("&")
        if len(parts) != 2:
            return None, None, None
        smrna_struct = parts[1]
        # Ensure the smRNA structure has sufficient length.
        seed_region = smrna_struct[1:8] if len(smrna_struct) >= 8 else smrna_struct[1:]
        central_region = smrna_struct[8:13] if len(smrna_struct) >= 13 else ""
        supplemental_region = smrna_struct[13:] if len(smrna_struct) > 13 else ""
        
        seed_bp = sum(1 for ch in seed_region if ch in "()")
        central_bp = sum(1 for ch in central_region if ch in "()")
        supplemental_bp = sum(1 for ch in supplemental_region if ch in "()")
        return seed_bp, central_bp, supplemental_bp
    except Exception as e:
        print(f"Error in annotate_bp_regions: {e}")
        return None, None, None

# -------------------------------
# Function to extract target interaction region (for genome annotation)
# -------------------------------
def extract_target_interaction_region(row, extension=2):
    """
    Extracts the region of the target sequence (extended_gene_mapped) that is predicted to interact with the smRNA.
    
    It examines the target portion of the Duplex Structure (before the "&"), finds the first and last indices
    where the structure is not ".", and then extends those boundaries by 'extension' nucleotides if possible.
    Returns the corresponding substring of the target sequence.
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

# Annotate base-pair counts in smRNA regions.
df_final[['Seed BP', 'Central BP', 'Supplemental BP']] = df_final.apply(
    lambda row: pd.Series(annotate_bp_regions(row["Duplex Structure"])), axis=1
)


def annotate_target_site_coordinates(row):
    """
    Annotates the target site coordinates so that the annotated region's length equals the length
    of the predicted target interaction region (TIR) extracted from the target side.
    
    The process is as follows:
      1. Use extract_target_interaction_region(row, extension=0) to obtain the TIR from the extended target sequence.
      2. Let L = len(TIR). Find the first occurrence of TIR in the extended target sequence.
      3. Since extended_gene_mapped was built as transcript[gene_start - 5 : gene_end + 5],
         index 0 corresponds to transcript coordinate (gene_start - 5).
         Thus, the preliminary genomic coordinates are:
              genomic_start = (gene_start - 5) + (index where TIR is found)
              genomic_end   = genomic_start + L - 1
      4. If these coordinates are off by 5 nt compared to your expectation, add 5 to both values.
    
    Returns a tuple (genomic_start, genomic_end) such that (genomic_end - genomic_start + 1) == L.
    """
    # Get the full target sequence used for duplex input.
    target_seq = row["extended_gene_mapped"]
    # Extract the predicted target interaction region (TIR) from the target side.
    tir = extract_target_interaction_region(row, extension=0)
    L = len(tir)
    if L == 0:
        # If nothing is predicted to pair, default to the entire gene_mapped region.
        return row["gene_start"], row["gene_end"]
    # Find the first occurrence of the TIR in the extended target sequence.
    idx = target_seq.find(tir)
    if idx == -1:
        idx = 0
        L = len(target_seq)
    # The extended target sequence was built using:
    #   start = gene_start - 5.
    # Thus, index 0 corresponds to transcript coordinate: gene_start - 5.
    # Preliminary mapping:
    offset = row["gene_start"] - 5
    preliminary_start = offset + idx
    preliminary_end = preliminary_start + L - 1
    # If your coordinates are off by 5 nt, add 5 to both.
    genomic_start = preliminary_start + 5
    genomic_end = preliminary_end + 5
    return genomic_start, genomic_end


# -------------------------------
# Apply the function to annotate target site coordinates.
# -------------------------------
df_final[['Target Site Start', 'Target Site End']] = df_final.apply(
    lambda row: pd.Series(annotate_target_site_coordinates(row)), axis=1
)

# Optionally, also compute the extracted target interaction region using the old function,
# so you can compare lengths.
df_final["Target Interaction Region"] = df_final.apply(extract_target_interaction_region, axis=1)

# For debugging, print the relevant columns:
print(df_final[["extended_gene_mapped", "Duplex Structure", "Target Interaction Region", 
                "Target Site Start", "Target Site End"]].head())

# Finally, save the updated DataFrame.
df_final.to_csv(output_file, sep="\t", index=False)
print(f"RNAduplex predictions complete! Results saved to {output_file}")
