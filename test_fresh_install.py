#!/usr/bin/env python3
"""
Test script for MolForge fresh installation.

This script tests basic functionality of MolForge in a fresh environment.
Run this after installing MolForge to verify the installation.

Usage:
    python test_fresh_install.py
"""

import sys
import pandas as pd
from pathlib import Path

def test_imports():
    """Test that all core modules can be imported."""
    print("=" * 70)
    print("Testing imports...")
    print("=" * 70)

    try:
        from molforge import (
            MolForge,
            ForgeParams,
            ConstructPipe,
            ChEMBLSourceParams,
            ChEMBLCuratorParams,
            CurateMolParams,
            TokenizeDataParams,
            CurateDistributionParams,
            GenerateConfsParams,
        )
        print("✓ All core imports successful")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False


def test_basic_curation():
    """Test basic molecule curation functionality."""
    print("\n" + "=" * 70)
    print("Testing basic molecule curation...")
    print("=" * 70)

    try:
        from molforge import MolForge, ForgeParams, CurateMolParams
        import tempfile

        # Create test data
        test_smiles = [
            'CCO',           # Ethanol
            'c1ccccc1',      # Benzene
            'CC(C)O',        # Isopropanol
            'CCN(CC)CC',     # Triethylamine
        ]

        test_df = pd.DataFrame({'SMILES': test_smiles})

        # Create temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            # Configure pipeline
            params = ForgeParams(
                steps=['curate'],
                output_root=tmpdir,
                write_checkpoints=False,
                curate_params=CurateMolParams(
                    SMILES_column='SMILES',
                    mol_steps=['sanitize', 'neutralize'],
                    smiles_steps=['canonical'],
                    dropna=True
                )
            )

            # Run pipeline
            forge = MolForge(params)
            result = forge.forge(test_df, input_id="test")

            # Verify results
            assert result is not None, "Result is None"
            assert isinstance(result, pd.DataFrame), "Result is not a DataFrame"
            assert len(result) > 0, "Result is empty"
            assert 'curated_smiles' in result.columns, "Missing curated_smiles column"

            print(f"✓ Curated {len(result)} molecules successfully")
            print(f"  Columns: {', '.join(result.columns)}")

        return True

    except Exception as e:
        print(f"✗ Curation test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_tokenization():
    """Test SMILES tokenization."""
    print("\n" + "=" * 70)
    print("Testing SMILES tokenization...")
    print("=" * 70)

    try:
        from molforge import MolForge, ForgeParams, CurateMolParams, TokenizeDataParams
        import tempfile

        test_smiles = ['CCO', 'c1ccccc1', 'CC(C)O']
        test_df = pd.DataFrame({'SMILES': test_smiles})

        with tempfile.TemporaryDirectory() as tmpdir:
            params = ForgeParams(
                steps=['curate', 'tokens'],
                output_root=tmpdir,
                write_checkpoints=False,
                curate_params=CurateMolParams(
                    SMILES_column='SMILES',
                    mol_steps=['sanitize'],
                    smiles_steps=['canonical'],
                    dropna=True
                ),
                tokens_params=TokenizeDataParams(
                    SMILES_column='curated_smiles',
                    dynamically_update_vocab=True
                )
            )

            forge = MolForge(params)
            result = forge.forge(test_df, input_id="test")

            assert 'tokens' in result.columns, "Missing tokens column"
            assert 'seqlen' in result.columns, "Missing seqlen column"

            print(f"✓ Tokenized {len(result)} molecules successfully")
            print(f"  Example tokens: {result['tokens'].iloc[0][:50]}...")

        return True

    except Exception as e:
        print(f"✗ Tokenization test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_conformer_generation():
    """Test 3D conformer generation."""
    print("\n" + "=" * 70)
    print("Testing conformer generation...")
    print("=" * 70)

    try:
        from molforge import MolForge, ForgeParams, CurateMolParams, GenerateConfsParams
        import tempfile

        # Use simple molecules for fast testing
        test_smiles = ['CCO', 'CC(C)O']
        test_df = pd.DataFrame({'SMILES': test_smiles})

        with tempfile.TemporaryDirectory() as tmpdir:
            params = ForgeParams(
                steps=['curate', 'confs'],
                output_root=tmpdir,
                write_checkpoints=False,
                curate_params=CurateMolParams(
                    SMILES_column='SMILES',
                    mol_steps=['sanitize'],
                    smiles_steps=['canonical'],
                    dropna=True
                ),
                confs_params=GenerateConfsParams(
                    backend='rdkit',
                    SMILES_column='curated_smiles',
                    max_confs=5,  # Small number for fast testing
                    rms_threshold=0.5,
                    use_uff=True,
                    max_iterations=100
                )
            )

            forge = MolForge(params)
            result = forge.forge(test_df, input_id="test")

            assert 'conformer_success' in result.columns, "Missing conformer_success column"
            assert 'n_conformers' in result.columns, "Missing n_conformers column"

            successful = result['conformer_success'].sum()
            print(f"✓ Generated conformers for {successful}/{len(result)} molecules")
            print(f"  Total conformers: {result['n_conformers'].sum()}")

        return True

    except Exception as e:
        print(f"✗ Conformer generation test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_construct_pipe():
    """Test ConstructPipe direct usage."""
    print("\n" + "=" * 70)
    print("Testing ConstructPipe...")
    print("=" * 70)

    try:
        from molforge import ConstructPipe, ForgeParams, CurateMolParams
        import tempfile

        test_smiles = ['CCO', 'c1ccccc1']
        test_df = pd.DataFrame({'SMILES': test_smiles})

        with tempfile.TemporaryDirectory() as tmpdir:
            params = ForgeParams(
                steps=['curate'],
                output_root=tmpdir,
                write_checkpoints=False,
                curate_params=CurateMolParams(
                    SMILES_column='SMILES',
                    mol_steps=['sanitize'],
                    smiles_steps=['canonical'],
                    dropna=True
                )
            )

            pipeline = ConstructPipe(params)
            result = pipeline(test_df, input_id="test_pipe")

            assert result is not None, "Pipeline result is None"
            assert len(result) > 0, "Pipeline result is empty"

            print(f"✓ ConstructPipe processed {len(result)} molecules successfully")

        return True

    except Exception as e:
        print(f"✗ ConstructPipe test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("MolForge Fresh Installation Test")
    print("=" * 70)

    tests = [
        ("Imports", test_imports),
        ("Basic Curation", test_basic_curation),
        ("Tokenization", test_tokenization),
        ("Conformer Generation", test_conformer_generation),
        ("ConstructPipe", test_construct_pipe),
    ]

    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n✗ Test '{name}' crashed: {e}")
            results.append((name, False))

    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed! MolForge is ready to use.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please check the output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
