def design_gibson_fragments(sequence, num_fragments=3, min_len=400, max_len=800, overlap=20):
    """
    Splits a DNA sequence into fragments suitable for Gibson Assembly.

    Parameters:
        sequence (str): DNA sequence (string of ACGT).
        num_fragments (int): Number of fragments to create (default = 3).
        min_len (int): Minimum fragment length.
        max_len (int): Maximum fragment length.
        overlap (int): Number of overlapping bases between fragments.

    Returns:
        list of dicts: Each dict has 'name', 'start', 'end', 'seq'
    """
    total_len = len(sequence)
    
    # Calculate ideal fragment length (excluding overlaps)
    total_overlap = (num_fragments - 1) * overlap
    available_length = total_len - total_overlap
    ideal_len = available_length // num_fragments
    
    # Ensure ideal length is within bounds
    if ideal_len < min_len:
        raise ValueError(f"Sequence too short for {num_fragments} fragments with {overlap}bp overlaps")
    
    # Clamp to max_len if needed
    frag_length = min(ideal_len, max_len)
    
    fragments = []
    start = 0
    
    for i in range(num_fragments):
        end = start + frag_length
        
        # For the last fragment, extend to the end of the sequence
        if i == num_fragments - 1:
            end = total_len
        
        fragment = {
            "name": f"Fragment {i + 1}",
            "start": start,
            "end": end,
            "seq": sequence[start:end]
        }
        fragments.append(fragment)
        
        # Next fragment starts with overlap
        start = end - overlap
    
    return fragments

# Example usage:
if __name__ == "__main__":
    def read_fasta(filepath):
        with open(filepath, "r") as f:
            lines = f.readlines()
        return "".join(line.strip() for line in lines if not line.startswith(">"))

    seq = read_fasta("pET28a_SHRT.fasta")
    frags = design_gibson_fragments(seq)
    for f in frags:
        print(f"{f['name']}: {f['start']}–{f['end']} ({len(f['seq'])} bp)")
        print(f"  Overlap preview: {f['seq'][:20]}...{f['seq'][-20:]}")
        print()
