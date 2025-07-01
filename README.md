![Banner](assets/github_banner.png)

# PCR Fragment Design

## Challenge 1: PCR Fragment Design for Gibson Assembly

Design 3 overlapping PCR fragments from a target DNA sequence. Each overlap must be 20 bp. Fragments must be between 400–800 bp in length.

### Requirements
- Input: `target_sequence.fasta` (single DNA string)
- Overlap length = 20 bp
- Min fragment length = 400 bp
- Max fragment length = 800 bp

### Output
- Coordinates + sequences of each fragment
- Check: overlaps match exactly

#### Example:
If the sequence is 2400 bp:
- Fragment 1: 0–~800 bp  
- Fragment 2: ~780–~1580 bp (shares 20 bp with Fragment 1)  
- Fragment 3: ~1560–2400 bp (shares 20 bp with Fragment 2)  

## Installation

1. Clone the repository:
```bash
git clone https://github.com/bmwoolf/PCR_Fragment_Design
cd PCR_Fragment_Design
```

2. Install dependencies:
```bash
pip install pytest
```

## Usage

### Using the Python Package

```python
from pcr_fragment_design.design_gibson_fragments import design_gibson_fragments

# Read your FASTA file
def read_fasta(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))

# Load your sequence
sequence = read_fasta("your_sequence.fasta")

# Design fragments
fragments = design_gibson_fragments(sequence)

# Print results
for fragment in fragments:
    print(f"{fragment['name']}: {fragment['start']}–{fragment['end']} ({len(fragment['seq'])} bp)")
    print(f"Sequence: {fragment['seq'][:50]}...")
    print()
```

### Custom Parameters

You can customize the fragment design parameters:

```python
fragments = design_gibson_fragments(
    sequence,
    num_fragments=4,    # Number of fragments (default: 3)
    min_len=300,        # Minimum fragment length (default: 400)
    max_len=600,        # Maximum fragment length (default: 800)
    overlap=15          # Overlap length (default: 20)
)
```

### Running the Demo

```bash
python scripts/run_tests.py
```

This will process the sample FASTA file and show detailed fragment information.

## Testing

Run all tests:
```bash
python -m pytest tests/ -v
```

Run a specific test:
```bash
python -m pytest tests/test_design_gibson_fragments.py::test_gibson_fragment_lengths -v
```

## Validation

The tests verify:
- Fragment lengths are within specified bounds
- Overlaps between adjacent fragments match exactly
- Fragment coordinates are correct
- Sequence extraction is accurate
- Error handling works properly

## Example Output

```
🧬 PCR Fragment Design Demo
==================================================
📁 Processing data/pET28a_SHRT.fasta...
📏 Sequence length: 2720 bp

🔬 Generated 3 fragments:
--------------------------------------------------
Fragment 1:
  📍 Coordinates: 0–800
  📏 Length: 800 bp
  🧬 Sequence (first 50bp): ATGACATCATTTAGTTGCCAGCCATGGTACAGTGAAAAGTTCTTCTCCTT...
  🔗 Overlap with Fragment 2: AGTTGCCAGCCATGGTACAG

✅ Overlap Verification:
  ✓ Fragment 1 ↔ Fragment 2: Overlap matches
  ✓ Fragment 2 ↔ Fragment 3: Overlap matches
```  