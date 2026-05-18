#!/usr/bin/env python3

from helpers import (
    reverse_complement,
    force_mismatch,
    parse_cigar_tuples,
    new_cigar,
    apply_n_mask,
    make_fastq_record,
)


def check_equal(test_name, observed, expected):
    if observed == expected:
        print(f"PASS: {test_name}")
    else:
        print(f"FAIL: {test_name}")
        print(f"  observed: {observed}")
        print(f"  expected: {expected}")
        raise AssertionError(test_name)


def check_raises(test_name, func, expected_message_part=None):
    try:
        func()
    except Exception as e:
        if expected_message_part is not None and expected_message_part not in str(e):
            print(f"FAIL: {test_name}")
            print(f"  raised message: {e}")
            print(f"  expected message to contain: {expected_message_part}")
            raise
        print(f"PASS: {test_name}")
    else:
        print(f"FAIL: {test_name}")
        print("  expected an error, but no error was raised")
        raise AssertionError(test_name)


def test_reverse_complement():
    check_equal(
        "reverse_complement basic",
        reverse_complement("ACGTN"),
        "NACGT",
    )

    check_equal(
        "reverse_complement lowercase",
        reverse_complement("acgtn"),
        "NACGT",
    )


def test_force_mismatch():
    check_equal("force_mismatch A", force_mismatch("A"), "G")
    check_equal("force_mismatch G", force_mismatch("G"), "A")
    check_equal("force_mismatch C", force_mismatch("C"), "T")
    check_equal("force_mismatch T", force_mismatch("T"), "C")
    check_equal("force_mismatch N", force_mismatch("N"), "N")
    check_equal("force_mismatch lowercase a", force_mismatch("a"), "G")


def test_parse_cigar_tuples():
    check_equal(
        "parse_cigar_tuples simple",
        parse_cigar_tuples("151M"),
        [(151, "M")],
    )

    check_equal(
        "parse_cigar_tuples complex",
        parse_cigar_tuples("10S50M2I20M1D71M"),
        [(10, "S"), (50, "M"), (2, "I"), (20, "M"), (1, "D"), (71, "M")],
    )

    check_raises(
        "parse_cigar_tuples missing",
        lambda: parse_cigar_tuples("*"),
        "Missing CIGAR",
    )

    check_raises(
        "parse_cigar_tuples malformed",
        lambda: parse_cigar_tuples("10MXYZ"),
        "Malformed CIGAR",
    )


def test_new_cigar_basic():
    # read: A C G T A C G T A A
    # ref:  A C G T T C G T A G
    #       = = = = X = = = = X
    check_equal(
        "new_cigar 10M with mismatches",
        new_cigar(
            cigar="10M",
            read_seq="ACGTACGTAA",
            ref_seq="ACGTTCGTAG",
        ),
        "4=1X4=1X",
    )

    check_equal(
        "new_cigar all matches",
        new_cigar(
            cigar="10M",
            read_seq="ACGTACGTAA",
            ref_seq="ACGTACGTAA",
        ),
        "10=",
    )

    check_equal(
        "new_cigar all mismatches",
        new_cigar(
            cigar="4M",
            read_seq="AAAA",
            ref_seq="TTTT",
        ),
        "4X",
    )


def test_new_cigar_no_m():
    # If there is no M, return unchanged.
    check_equal(
        "new_cigar no M returns unchanged",
        new_cigar(
            cigar="5=1X4=",
            read_seq="ACGTACGTAA",
            ref_seq="ACGTTCGTAA",
        ),
        "5=1X4=",
    )

    check_equal(
        "new_cigar no M with softclip returns unchanged",
        new_cigar(
            cigar="10S141=",
            read_seq="A" * 151,
            ref_seq="A" * 141,
        ),
        "10S141=",
    )


def test_new_cigar_with_indels_and_clipping():
    # CIGAR: 5M2I5M1D5M
    #
    # Query-consuming:
    #   5M + 2I + 5M + 5M = 17 query bases
    #
    # Reference-consuming:
    #   5M + 5M + 1D + 5M = 16 reference bases
    #
    # First 5M:
    #   read AAAAA vs ref AAAAA => 5=
    #
    # 2I:
    #   preserved as 2I
    #
    # Next 5M:
    #   read GGGGA vs ref GGGGC => 4=1X
    #
    # 1D:
    #   preserved as 1D
    #
    # Final 5M:
    #   read TTTTT vs ref TTTTT => 5=
    check_equal(
        "new_cigar with insertion and deletion",
        new_cigar(
            cigar="5M2I5M1D5M",
            read_seq="AAAAACCGGGGATTTTT",
            ref_seq="AAAAAGGGGCGTTTTT",
        ),
        "5=2I4=1X1D5=",
    )

    # Soft clip consumes query but not reference.
    check_equal(
        "new_cigar with softclip",
        new_cigar(
            cigar="3S5M",
            read_seq="NNNACGTA",
            ref_seq="ACGTT",
        ),
        "3S4=1X",
    )


def test_new_cigar_errors():
    check_raises(
        "new_cigar None cigar",
        lambda: new_cigar(None, "ACGT", "ACGT"),
        "Cannot create new CIGAR",
    )

    check_raises(
        "new_cigar star cigar",
        lambda: new_cigar("*", "ACGT", "ACGT"),
        "Cannot create new CIGAR",
    )

    check_raises(
        "new_cigar read too short",
        lambda: new_cigar("10M", "ACGT", "ACGTACGTAC"),
        "read sequence is shorter",
    )

    check_raises(
        "new_cigar ref too short",
        lambda: new_cigar("10M", "ACGTACGTAC", "ACGT"),
        "reference sequence is shorter",
    )


def test_apply_n_mask():
    check_equal(
        "apply_n_mask basic",
        apply_n_mask(
            new_seq="ACGTACGTAA",
            original_seq="AAGNACNTAA",
        ),
        "ACGNACNTAA",
    )

    check_equal(
        "apply_n_mask lowercase original",
        apply_n_mask(
            new_seq="ACGTACGTAA",
            original_seq="aagnacntaa",
        ),
        "ACGNACNTAA",
    )

    check_raises(
        "apply_n_mask length mismatch",
        lambda: apply_n_mask("ACGT", "ACGTA"),
        "sequence lengths differ",
    )


def test_make_fastq_record():
    expected = "@read1/1\nACGT\n+\nIIII\n"

    check_equal(
        "make_fastq_record default Q40",
        make_fastq_record("read1/1", "ACGT"),
        expected,
    )

    expected_custom = "@read1/1\nACGT\n+\n!!!!\n"

    check_equal(
        "make_fastq_record custom quality",
        make_fastq_record("read1/1", "ACGT", phred_char="!"),
        expected_custom,
    )


def main():
    print("Testing helpers.py functions...\n")

    test_reverse_complement()
    test_force_mismatch()
    test_parse_cigar_tuples()
    test_new_cigar_basic()
    test_new_cigar_no_m()
    test_new_cigar_with_indels_and_clipping()
    test_new_cigar_errors()
    test_apply_n_mask()
    test_make_fastq_record()

    print("\nAll helper tests passed.")


if __name__ == "__main__":
    main()