#!/usr/bin/env python3
"""
PCR Fragment Design Test Runner
Runs all tests and demonstrates the functionality
"""

import sys
import os
from pcr_fragment_design.design_gibson_fragments import design_gibson_fragments

def read_fasta(filepath):
    """Helper function to read FASTA files"""
    with open(filepath, "r") as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))

def run_demo():
    """Run a demonstration of the PCR fragment design"""
    print("🧬 PCR Fragment Design Demo")
    print("=" * 50)
    
    # Test with the real FASTA file
    fasta_path = os.path.join("data", "pET28a_SHRT.fasta")
    if os.path.exists(fasta_path):
        print(f"📁 Processing {fasta_path}...")
        seq = read_fasta(fasta_path)
        print(f"📏 Sequence length: {len(seq)} bp")
        
        frags = design_gibson_fragments(seq)
        
        print(f"\n🔬 Generated {len(frags)} fragments:")
        print("-" * 50)
        
        for i, frag in enumerate(frags):
            print(f"{frag['name']}:")
            print(f"  📍 Coordinates: {frag['start']}–{frag['end']}")
            print(f"  📏 Length: {len(frag['seq'])} bp")
            print(f"  🧬 Sequence (first 50bp): {frag['seq'][:50]}...")
            print(f"  🧬 Sequence (last 50bp): ...{frag['seq'][-50:]}")
            
            # Show overlap with next fragment
            if i < len(frags) - 1:
                next_frag = frags[i + 1]
                overlap = frag['seq'][-20:]
                print(f"  🔗 Overlap with {next_frag['name']}: {overlap}")
            print()
        
        # Verify overlaps
        print("✅ Overlap Verification:")
        for i in range(len(frags) - 1):
            current_end = frags[i]['seq'][-20:]
            next_start = frags[i+1]['seq'][:20]
            if current_end == next_start:
                print(f"  ✓ {frags[i]['name']} ↔ {frags[i+1]['name']}: Overlap matches")
            else:
                print(f"  ✗ {frags[i]['name']} ↔ {frags[i+1]['name']}: Overlap mismatch!")
                print(f"    Expected: {current_end}")
                print(f"    Got: {next_start}")
    
    else:
        print("❌ FASTA file not found. Creating test sequence...")
        # Create a test sequence
        test_seq = "ATCG" * 600  # 2400 bp
        frags = design_gibson_fragments(test_seq)
        
        print(f"📏 Test sequence length: {len(test_seq)} bp")
        print(f"🔬 Generated {len(frags)} fragments:")
        
        for frag in frags:
            print(f"  {frag['name']}: {frag['start']}–{frag['end']} ({len(frag['seq'])} bp)")

def run_basic_tests():
    """Run basic functionality tests"""
    print("\n🧪 Running Basic Tests")
    print("=" * 30)
    
    # Test 1: Basic functionality
    print("Test 1: Basic fragment generation...")
    seq = "ATCG" * 600  # 2400 bp
    frags = design_gibson_fragments(seq)
    
    assert len(frags) == 3, f"Expected 3 fragments, got {len(frags)}"
    print("  ✓ Correct number of fragments")
    
    # Test 2: Fragment lengths
    for i, frag in enumerate(frags[:-1]):  # Last fragment can be longer
        assert 400 <= len(frag['seq']) <= 800, f"Fragment {i+1} length {len(frag['seq'])} outside range"
    print("  ✓ Fragment lengths within range")
    
    # Test 3: Overlaps
    for i in range(len(frags) - 1):
        assert frags[i]['seq'][-20:] == frags[i+1]['seq'][:20], f"Overlap mismatch between fragments {i+1} and {i+2}"
    print("  ✓ Overlaps match exactly")
    
    # Test 4: Coordinates
    for i in range(len(frags) - 1):
        overlap = frags[i]['end'] - frags[i+1]['start']
        assert overlap == 20, f"Expected 20bp overlap, got {overlap}bp"
    print("  ✓ Coordinates correct")
    
    print("✅ All basic tests passed!")

if __name__ == "__main__":
    try:
        run_demo()
        run_basic_tests()
        print("\n🎉 All tests completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1) 