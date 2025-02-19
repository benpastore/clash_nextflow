import pandas as pd
import subprocess
from multiprocessing import Pool
import re
from Bio import SeqIO  # For parsing FASTA files
import sys
import multiprocessing
from bear_encoder.BEAR_encoder.bear import annotate_bp_bear_encoder

# -------------------------------
# (1) Load FASTA sequences into dictionaries
# -------------------------------
def load_fasta_sequences(fasta_file):
    """Parses a FASTA file and returns a dictionary {sequence_id: sequence}."""
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = str(record.seq).upper()  # Convert to uppercase
    return fasta_dict

# -------------------------------
# Functions to process sequences
# -------------------------------
def extend_gene_mapped(row, transcript_sequences):
    """Extend the target sequence (gene_mapped) by ±5 nt using gene_start and gene_end."""
    gene_id = row["gene_id"]
    if gene_id in transcript_sequences:
        full_seq = transcript_sequences[gene_id]
        start = max(0, row["gene_start"] - 5)
        end = min(len(full_seq), row["gene_end"] + 5)
        return full_seq[start:end]
    return row["gene_mapped"]

def fix_smrna_sequence(row, smrna_sequences):
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

# -------------------------------
# Bulge detection function 
# -------------------------------
def detect_bulges(dot_bracket):
    """
    Detects "true bulges" in a dot-bracket structure and ensures single nucleotide mismatches are not called as bulges.

    A region is classified as a bulge if:
    - It is at least **1 nt long**.
    - It is **flanked by base pairs on at least one strand**.
    - The **complementary strand has at least one paired base adjacent to the bulge**.

    Returns a tuple: (target_bulges, smrna_bulges)
    """
    if "&" not in dot_bracket:
        return "Invalid Structure", "Invalid Structure"

    try:
        split_index = dot_bracket.index("&")
        target_structure = dot_bracket[:split_index]  # Target RNA structure
        smrna_structure = dot_bracket[split_index+1:]  # smRNA structure

        def find_true_bulges(structure, complementary_structure, offset=0):
            bulge_positions = []
            L = len(structure)
            L_comp = len(complementary_structure)

            for match in re.finditer(r"\.+", structure):  # Find unpaired regions
                start, end = match.start(), match.end() - 1  # Get bulge indices

                # Ensure bulge is at least 1 nt long
                region_length = end - start + 1
                if region_length < 1:
                    continue  # Ignore single mismatches

                # Ensure bulge is NOT at sequence boundaries
                if start == 0 or end == L - 1:
                    continue  # Ignore boundary mismatches

                # **Check flanking base pairs on the same strand**
                left_paired = (start > 0 and structure[start - 1] in "()")  
                right_paired = (end + 1 < L and structure[end + 1] in "()")  

                # **Check flanking base pairs on the complementary strand**
                left_comp_paired = (start > 0 and start - 1 < L_comp and complementary_structure[start - 1] in "()")
                right_comp_paired = (end + 1 < L_comp and complementary_structure[end + 1] in "()")

                # A valid bulge must be flanked by base pairs on at least ONE strand
                if not (left_paired or right_paired or left_comp_paired or right_comp_paired):
                    continue  # Ignore random unpaired segments

                # **Ensure at least one paired base in the complementary strand**
                effective_start = max(0, start)
                effective_end = min(end, L_comp - 1)

                paired_bases = sum(1 for i in range(effective_start, effective_end + 1)
                                   if i < len(complementary_structure) and complementary_structure[i] in "()")

                if paired_bases > 0:
                    bulge_positions.append(f"Bulge at {start+offset}-{end+offset}")

            return ", ".join(bulge_positions) if bulge_positions else "No Bulges"

        # Detect bulges separately for target and smRNA
        target_bulges = find_true_bulges(target_structure, smrna_structure)
        smrna_bulges = find_true_bulges(smrna_structure, target_structure, offset=0)

        return target_bulges, smrna_bulges

    except Exception as e:
        return f"Error processing structure: {e}", "Error"

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
    
    target_seq, smrna_seq = [s.strip() for s in seq_pair.split("&")]
    duplex_input = f"{target_seq}\n{smrna_seq}\n"
    result = subprocess.run(
        ["RNAduplex", "-T", "20"],
        input=duplex_input,
        capture_output=True,
        text=True
    )
    output = result.stdout.strip()

    line = output  
    entries = [ i for i in line.strip().split(" ") if not i == '' and not i == ':' ]
    structure = entries[0]
    target_structure = structure.split("&")[0]
    smRNA_structure = structure.split("&")[1]
    target_interaction_residues = entries[1]
    smrna_interaction_residues = entries[2]
    binding_energy = entries[3].replace('(', "").replace(")", "")

    target_seq_len = len(target_seq)
    target_interaction_start = int( target_interaction_residues.split(",")[0] )
    target_interaction_end = int( target_interaction_residues.split(",")[1] )
    target_structure_fixed = "."*(target_interaction_start - 1) + target_structure + "."*(target_seq_len - target_interaction_end)
    
    # get an estimate of the sequence targeted by piRNA
    ts_start = target_interaction_end - 22 if target_interaction_end - 22 > 0 else 0
    ts_end = target_interaction_end 
    ts_seq = target_seq[ts_start:ts_end]

    smrna_seq_len = len(smrna_seq)
    smrna_interaction_start = int( smrna_interaction_residues.split(",")[0] )
    smrna_interaction_end = int( smrna_interaction_residues.split(",")[1] )
    smrna_structure_fixed = "."*(smrna_interaction_start - 1) + smRNA_structure + "."*(smrna_seq_len - smrna_interaction_end)
    
    structure_fixed = target_structure_fixed + "&" + smrna_structure_fixed

    bear = annotate_bp_bear_encoder(seq_pair, structure_fixed)
    return seq_pair, binding_energy, structure_fixed, ts_seq, bear

# -------------------------------
# Tail base-pairing analysis function
# -------------------------------
def is_tail_base_paired(row):

    tail = row['tail']
    if tail == "*" : 
        return "none"
    
    smrna_bp = row['structure_fixed'].split("&")[1]
    tail_len = len(tail)
    is_bp = True if ")" in smrna_bp[-1*tail_len:] else False

    return is_bp

# -------------------------------
# New Function: Annotate Base-Pair Counts in smRNA Regions
# -------------------------------
def annotate_bp_regions(duplex_structure):

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
    

def main() : 

    # -------------------------------
    # File paths
    # -------------------------------

    #input_file = sys.argv[1]
    #transcript_fasta = sys.argv[2]
    #smrna_fasta = sys.argv[3]
    #output_file = sys.argv[4]

    input_file = "/fs/ess/PCON0160/ben/pipelines/clash_nextflow/test.tsv"  # TSV file with gene coordinates
    transcript_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/ce_ws279.linc.pseudo.pc.repbase.fa"  # Transcript sequences
    smrna_fasta = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNA_reference.15ntextension.fa"  # smRNA sequences
    output_file = "RNAduplex_results.tsv"

    transcript_sequences = load_fasta_sequences(transcript_fasta)
    smrna_sequences = load_fasta_sequences(smrna_fasta)

    df = pd.read_csv(input_file, sep="\t")
    df = df.head(100)

    # Process sequences and convert T -> U.
    df["extended_gene_mapped"] = df.apply(lambda x: extend_gene_mapped(x, transcript_sequences), axis=1)
    df["corrected_smRNA_seq"] = df.apply(lambda x: fix_smrna_sequence(x, smrna_sequences), axis=1)

    df["extended_gene_mapped"] = df["extended_gene_mapped"].str.replace("T", "U")
    df["corrected_smRNA_seq"] = df["corrected_smRNA_seq"].str.replace("T", "U")

    # Create RNAduplex input in the form: "extended_gene_mapped&corrected_smRNA_seq"
    df["duplex_input"] = df["extended_gene_mapped"] + "&" + df["corrected_smRNA_seq"]

    # Remove duplicate duplex inputs.
    df_unique = df.drop_duplicates(subset=["duplex_input"]).copy()

    NUM_CORES = multiprocessing.cpu_count()-2  # Adjust as needed
    with Pool(NUM_CORES) as p:
        results = p.map(run_rnaduplex, df_unique["duplex_input"].tolist())

    # Create a DataFrame from RNAduplex results.
    df_results = pd.DataFrame(results, columns=["duplex_input", "binding_energy", "structure_fixed", "target_site_seq", "bear"])


    # Merge RNAduplex results with the original DataFrame.
    df_final = df.merge(df_results, on="duplex_input", how="left")

    # Add Tail Base-Paired Analysis.
    df_final["tail_bp"] = df_final.apply(is_tail_base_paired, axis=1)
    df_final.to_csv(output_file, sep="\t", index=False)


    # Annotate base-pair counts in smRNA regions.
    df_final[['seed_bp', 'central_bp', 'supplemental_bp']] = df_final.apply(
        lambda row: pd.Series(annotate_bp_regions(row["structure_fixed"])), axis=1
    )

    # Finally, save the updated DataFrame.
    df_final.to_csv(output_file, sep="\t", index=False)
    print(f"RNAduplex predictions complete! Results saved to {output_file}")


if __name__ == "__main__" : 

    main()









