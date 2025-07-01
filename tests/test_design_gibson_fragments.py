import pytest
from pcr_fragment_design.design_gibson_fragments import design_gibson_fragments

def read_fasta(filepath):
    """Helper function to read FASTA files"""
    with open(filepath, "r") as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))

def test_gibson_fragment_lengths():
    """Test basic fragment length requirements"""
    seq = "A" * 2400  # mock DNA
    frags = design_gibson_fragments(seq)

    assert len(frags) == 3
    for i in range(2):
        # Check 20 bp overlaps
        assert frags[i]['seq'][-20:] == frags[i+1]['seq'][:20]
        # Check length limits
        assert 400 <= len(frags[i]['seq']) <= 800

def test_overlap_verification():
    """Test that overlaps between fragments match exactly"""
    seq = "ATCG" * 600  # 2400 bp sequence with repeating pattern
    frags = design_gibson_fragments(seq)
    
    # Check overlaps between adjacent fragments
    for i in range(len(frags) - 1):
        current_frag_end = frags[i]['seq'][-20:]
        next_frag_start = frags[i+1]['seq'][:20]
        assert current_frag_end == next_frag_start, f"Overlap mismatch between fragments {i+1} and {i+2}"

def test_fragment_coordinates():
    """Test that fragment coordinates are correct and non-overlapping (except for overlaps)"""
    seq = "ATCG" * 600
    frags = design_gibson_fragments(seq)
    
    for i in range(len(frags) - 1):
        # Check that fragments don't overlap beyond the specified overlap length
        current_end = frags[i]['end']
        next_start = frags[i+1]['start']
        actual_overlap = current_end - next_start
        assert actual_overlap == 20, f"Expected 20bp overlap, got {actual_overlap}bp"

def test_fragment_sequence_extraction():
    """Test that fragment sequences are correctly extracted from the original sequence"""
    seq = "ATCG" * 600
    frags = design_gibson_fragments(seq)
    
    for frag in frags:
        extracted_seq = seq[frag['start']:frag['end']]
        assert frag['seq'] == extracted_seq, f"Sequence mismatch for {frag['name']}"

def test_real_fasta_file():
    """Test with the actual FASTA file"""
    try:
        seq = read_fasta("data/pET28a_SHRT.fasta")
        frags = design_gibson_fragments(seq)
        
        # Basic checks
        assert len(frags) == 3
        assert all(400 <= len(f['seq']) <= 800 for f in frags[:-1])  # Last fragment can be longer
        
        # Check overlaps
        for i in range(len(frags) - 1):
            assert frags[i]['seq'][-20:] == frags[i+1]['seq'][:20]
            
        print(f"Successfully processed {len(seq)}bp sequence into {len(frags)} fragments")
        for f in frags:
            print(f"{f['name']}: {f['start']}–{f['end']} ({len(f['seq'])} bp)")
            
    except FileNotFoundError:
        pytest.skip("FASTA file not found")

def test_edge_case_minimum_sequence():
    """Test that a ValueError is raised for minimum sequence length"""
    min_seq = "ATCG" * 290  # 1160 bp
    with pytest.raises(ValueError):
        design_gibson_fragments(min_seq)

def test_edge_case_maximum_fragment_length():
    """Test when fragments should be at maximum length"""
    # Sequence long enough to hit max fragment length
    seq = "ATCG" * 1000  # 4000 bp
    frags = design_gibson_fragments(seq)
    
    # First two fragments should be at max length (800bp)
    assert len(frags[0]['seq']) == 800
    assert len(frags[1]['seq']) == 800
    # Last fragment can be longer to cover remaining sequence

def test_custom_parameters():
    """Test with custom fragment parameters"""
    seq = "ATCG" * 400  # 1600 bp
    frags = design_gibson_fragments(seq, num_fragments=4, min_len=300, max_len=500, overlap=15)
    
    assert len(frags) == 4
    assert all(300 <= len(f['seq']) <= 500 for f in frags[:-1])
    
    # Check custom overlap length
    for i in range(len(frags) - 1):
        assert frags[i]['seq'][-15:] == frags[i+1]['seq'][:15]

def test_sequence_too_short():
    """Test error handling for sequences too short"""
    short_seq = "ATCG" * 50  # 200 bp - too short for 3 fragments with 20bp overlaps
    with pytest.raises(ValueError):
        design_gibson_fragments(short_seq)

def test_dna_sequence_validation():
    """Test that a ValueError is raised for too-short DNA sequences"""
    valid_dna = "ATCGATCGATCG"
    with pytest.raises(ValueError):
        design_gibson_fragments(valid_dna, min_len=3, max_len=6, overlap=2)

if __name__ == "__main__":
    # Run all tests
    pytest.main([__file__, "-v"])
