# PCR_Fragment_Design

## Challenge 1: PCR Fragment Design for Gibson Assembly

Design 3 overlapping PCR fragments from a target DNA sequence. Each overlap must be 20 bp. Fragments must be between 400–800 bp in length.

Input:
- `target_sequence.fasta` (single DNA string)
- Overlap length = 20 bp
- Min fragment length = 400 bp
- Max fragment length = 800 bp

Output:
- Coordinates + sequences of each fragment
- Check: overlaps match exactly

#### Example:
If the sequence is 2400 bp:
- Fragment 1: 0–~800 bp  
- Fragment 2: ~780–~1580 bp (shares 20 bp with Fragment 1)  
- Fragment 3: ~1560–2400 bp (shares 20 bp with Fragment 2)  