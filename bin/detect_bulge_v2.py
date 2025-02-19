def detect_bulges_rnaduplex_interaction(full_seq_target,
                                        full_seq_smrna,
                                        full_dot_bracket):
    """
    Detect bulges from an RNA–RNA interaction dot-bracket (RNAduplex style).
    
    :param full_seq_target: String of the target RNA (same length as left side of dot-bracket).
    :param full_seq_smrna:  String of the smRNA (same length as right side of dot-bracket).
    :param full_dot_bracket: Dot-bracket structure for the interaction, including '&'.
                             Example: "...........((((((((((((((((....&)))))).)))))).))))."
    :return: {
       "pairs":       List of (t_i, s_i) pairs,
       "bulges_t":    List of (start, end) bulges on the target (indices in target_seq),
       "bulges_s":    List of (start, end) bulges on the smRNA (indices in smRNA),
       "maybe_reversed": Boolean indicating if smRNA was reversed
    }
    """
    # 1) Split at '&' to find the boundary.
    #    Left portion => target structure, Right portion => smRNA structure
    #    BUT for an inter-molecular fold, parentheses often cross from left to right.
    #    So we must parse them as ONE string, ignoring '&' in the stack logic.
    
    if '&' not in full_dot_bracket:
        raise ValueError("No '&' found in dot-bracket for two-strand interaction.")
    
    mid = full_dot_bracket.index('&')
    left_part  = full_dot_bracket[:mid]
    right_part = full_dot_bracket[mid+1:]
    
    # Check if the right part might be reversed (common if it begins with ')')
    # This is a heuristic; adjust if you have a more reliable check.
    maybe_reversed = False
    if right_part.strip().startswith(')'):
        # Reverse the smRNA sequence and its bracket portion
        full_seq_smrna = full_seq_smrna[::-1]
        right_part     = right_part[::-1]
        maybe_reversed = True
    
    # Reconstruct the full bracket string without '&' so we can parse pairs in one pass
    combined_bracket = left_part + right_part
    
    # Sanity check lengths:
    # Left bracket length should match len(full_seq_target).
    # Right bracket length should match len(full_seq_smrna).
    if len(left_part)  != len(full_seq_target):
        raise ValueError("Left dot-bracket length does not match target sequence length.")
    if len(right_part) != len(full_seq_smrna):
        raise ValueError("Right dot-bracket length does not match smRNA sequence length.")
    
    # 2) Build the list of base pairs using a stack approach
    stack = []
    pairs = []  # will store tuples (index_in_combined, index_in_combined)
    
    for i, char in enumerate(combined_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                # Mismatched parenthesis
                continue
            opening_index = stack.pop()
            pairs.append((opening_index, i))
        # If it's '.', just do nothing
    
    # Now we have pairs in terms of combined indices [0..(len(combined_bracket)-1)].
    # Next, we must figure out which pairs are inter-molecular:
    #   - target indices: 0..(len(left_part)-1)
    #   - smRNA indices:  0..(len(right_part)-1)
    # but offset by mid+1 in the original notation. We'll map them carefully.

    # Create a list of (t_i, s_i) for each pair that crosses from left to right.
    # If a pair is (x, y) with x < len(left_part) and y >= len(left_part),
    # that means x is on target side, y is on smRNA side.
    # We'll map:
    #   target_index = x
    #   smRNA_index  = y - len(left_part)
    # (Similarly if x > left_part => smRNA side, but typically for an RNA–RNA duplex
    #  the '(' appear on the left side and the ')' appear on the right side.)
    
    inter_pairs = []
    L = len(left_part)
    for (open_i, close_i) in pairs:
        # Sort so open_i < close_i for clarity
        i1, i2 = sorted([open_i, close_i])
        if i1 < L <= i2:
            # i1 is in the target, i2 is in the smRNA
            t_idx = i1
            s_idx = i2 - L
            inter_pairs.append((t_idx, s_idx))
        # If both i1, i2 < L, that's intramolecular pairing in the target
        # If both i1, i2 >= L, that's intramolecular pairing in the smRNA
        # We ignore intramolecular pairs for bulge detection between strands.
    
    # 3) Sort inter_pairs by the target index ascending
    inter_pairs.sort(key=lambda x: x[0])
    
    # 4) Detect bulges by scanning consecutive pairs
    bulges_target = []
    bulges_smrna  = []
    
    for i in range(len(inter_pairs) - 1):
        (t1, s1) = inter_pairs[i]
        (t2, s2) = inter_pairs[i+1]
        
        dt = t2 - t1
        ds = s2 - s1
        
        # If dt > 1 and ds == 1 => bulge on target side
        if dt > 1 and ds == 1:
            # Bulge is the region t1+1..t2-1 on the target
            bulges_target.append((t1+1, t2-1))
        
        # If ds > 1 and dt == 1 => bulge on smRNA side
        if ds > 1 and dt == 1:
            # Bulge is the region s1+1..s2-1 on the smRNA
            bulges_smrna.append((s1+1, s2-1))
        
        # If both dt>1 and ds>1 => interior loop (bulge on both sides)
        # If you want to record interior loops, you can do:
        # if dt > 1 and ds > 1:
        #    bulges_target.append((t1+1, t2-1))
        #    bulges_smrna.append((s1+1, s2-1))
    
    return {
        "pairs": inter_pairs,              # List of (t_i, s_i)
        "bulges_t": bulges_target,         # Bulges on the target side
        "bulges_s": bulges_smrna,          # Bulges on the smRNA side
        "maybe_reversed": maybe_reversed,  # Whether we reversed the smRNA
    }

# -------------------- EXAMPLE USAGE --------------------

if __name__ == "__main__":
    # From your example:
    # Full target (31 nt)  : CUCUUCCACUCGUUUAUUCUUAUGUGAGCAA
    # Full smRNA (19 nt)   : UCACAUCAAGAAUUGAACA
    # Dot-bracket (31 + 19 = 50 chars + '&'):
    #   "...........((((((((((((((((....&)))))).)))))).))))."
    #
    # Where the left part has length 31, the right part has length 19.

    target_seq = "CGUUUAUUCUUAUGUGAG"
    smrna_seq  = "ACAAGUUAAGAACUACACU"
    db_struct  = ".((((((((((((((((.&)))))).)))))).))))."

    result = detect_bulges_rnaduplex_interaction(target_seq, smrna_seq, db_struct)

    print("Inter-molecular base pairs (target_idx, smRNA_idx):")
    for p in result["pairs"]:
        print("  ", p)

    print("\nBulges on Target:", result["bulges_t"])
    print("Bulges on smRNA:", result["bulges_s"])
    print("Was smRNA reversed?", result["maybe_reversed"])
