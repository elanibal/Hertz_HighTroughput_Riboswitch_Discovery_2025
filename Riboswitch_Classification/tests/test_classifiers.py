"""Offline unit tests for the riboswitch classifiers (no network).

Run:  python -m pytest tests/ -q      (or)   python tests/test_classifiers.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from classify.terminator import find_terminator
from classify.sd_detector import detect_sd
from classify.direction import predict_direction
from classify.classifier import classify_sequence
from classify.microorf import find_microorfs
from classify import schema

# --- validated synthetic building blocks ---------------------------------------
STRONG_HAIRPIN = "GGGGGCGCGAAAACGCGCCCCC"      # stem + AAAA loop
POLYU = "TTTTTTT"
EP_TERMINATOR = "ATGCATGCATGCATGC" + STRONG_HAIRPIN + POLYU + "ACGTACGT"

ANTI_SD_TAIL = "TTTTTTCCTCCTTTTTTT"           # CCTCCT pairs with AGGAGG
EP_TRANSLATIONAL = "AGGAGGA" + "AAAAAAAA" + "ATG" + "AAACGT"  # SD@0, ATG@15
EP_NEITHER = "ACACACACACACAGCTGACGATCGATCGATCGA"


# --- terminator ----------------------------------------------------------------
def test_terminator_found():
    h = find_terminator(EP_TERMINATOR)
    assert h.found is True
    assert float(h.dG) < -8.0
    assert h.polyU.startswith("TTTTT")
    assert h.coords()

def test_terminator_polyU_without_hairpin_is_not_terminator():
    h = find_terminator("ACACACACACACACACACAC" + POLYU + "ACAC")
    assert h.found is False
    assert "no stable upstream hairpin" in h.reason

def test_terminator_absent_without_polyU():
    h = find_terminator("ACGACGACGACGACGACGACG")
    assert h.found is False
    assert h.reason == "no poly-U tract"


# --- SD detection / sequestration ----------------------------------------------
def test_sd_sequestered():
    sd = detect_sd(EP_TRANSLATIONAL, nts_to_start_codon=15, aptamer_seq=ANTI_SD_TAIL)
    assert sd.found is True
    assert sd.motif == "AGGAGG"
    assert sd.sequestered is True
    assert float(sd.paired_fraction) >= 0.5

def test_sd_present_but_accessible():
    sd = detect_sd(EP_TRANSLATIONAL, nts_to_start_codon=15, aptamer_seq="A" * 30)
    assert sd.found is True
    assert sd.sequestered is False           # no complementary tail -> not sequestered

def test_sd_absent():
    sd = detect_sd("CCCCCCCCCCCCATG", nts_to_start_codon=12, aptamer_seq="AAAA")
    assert sd.found is False


# --- direction rubric ----------------------------------------------------------
def test_direction_off_biosynthesis():
    d = predict_direction("thiamine biosynthesis protein ThiC", "thiC", "RF00059")
    assert d.direction == "OFF" and d.basis == "annotation"

def test_direction_on_transport():
    d = predict_direction("thiamine ABC transporter permease", "thiT", "RF00059")
    assert d.direction == "ON"

def test_direction_unknown_hypothetical():
    d = predict_direction("hypothetical protein", "", "RF00059")
    assert d.direction == "unknown" and d.basis == "none"

def test_direction_generic_transport():
    d = predict_direction("MFS transporter", "", "RF99999")
    assert d.direction == "ON"


# --- end-to-end row + schema integrity -----------------------------------------
def test_row_transcriptional_off():
    row = classify_sequence(
        rfam_family="RF00059", aptamer_id="TEST/1-50", accession="TEST",
        genome_coords="1-50(+)", strand="+", aptamer_seq="A" * 50,
        ep_seq=EP_TERMINATOR, downstream_gene_product="thiamine biosynthesis protein ThiC",
        downstream_gene="thiC", nts_to_start_codon=120,
    )
    assert row["mechanism_call"] == "transcriptional"
    assert row["predicted_direction"] == "OFF"
    assert row["has_intrinsic_terminator"] is True
    assert schema.validate_row(row) == []

def test_row_translational():
    # aptamer supplies the anti-SD tail that sequesters the SD
    row = classify_sequence(
        rfam_family="RF00059", aptamer_id="TEST/1-18", accession="TEST",
        genome_coords="1-18(+)", strand="+", aptamer_seq=ANTI_SD_TAIL,
        ep_seq=EP_TRANSLATIONAL, downstream_gene_product="thiamine transporter",
        downstream_gene="thiT", nts_to_start_codon=15,
    )
    assert row["mechanism_call"] == "translational"
    assert row["has_SD"] is True and row["anti_SD_overlap"] is True
    assert schema.validate_row(row) == []

def test_row_other_kept_not_dropped():
    row = classify_sequence(
        rfam_family="RF00059", aptamer_id="TEST/1-33", accession="TEST",
        genome_coords="1-33(+)", strand="+", aptamer_seq="A" * 33,
        ep_seq=EP_NEITHER, downstream_gene_product="", nts_to_start_codon="",
    )
    assert row["mechanism_call"] == "other"
    assert row["notes"]                       # a reason is recorded (no silent drop)
    assert schema.validate_row(row) == []

def test_all_columns_present_and_ordered():
    row = classify_sequence(
        rfam_family="RF00059", aptamer_id="x", accession="x", genome_coords="x",
        strand="+", aptamer_seq="A" * 20, ep_seq=EP_TERMINATOR,
    )
    assert list(row.keys()) == schema.COLUMNS  # canonical order preserved

def test_microorf_finds_small_orf():
    # ATG + 4 codons + TAA  (aa_len 5), embedded in flanks
    leader = "CCCC" + "ATG" + "AAAGGGCATTGC" + "TAA" + "GGGG"
    res = find_microorfs(leader, min_aa=2, max_aa=50)
    assert res.count >= 1
    assert any(o.start_codon == "ATG" and o.aa_len == 5 for o in res.orfs)

def test_microorf_respects_length_bounds():
    leader = "ATG" + "TAA"          # zero-codon ORF -> below min_aa
    assert find_microorfs(leader, min_aa=2).count == 0

def test_microorf_none_when_no_orf():
    assert find_microorfs("AAAAAAAAAAAAAAAA").count == 0

def test_row_has_microorf_and_truncation_fields():
    row = classify_sequence(
        rfam_family="RF00059", aptamer_id="x", accession="x", genome_coords="x",
        strand="+", aptamer_seq="A" * 20, ep_seq=EP_TERMINATOR,
        cds_extends_beyond_window=True,
    )
    assert "microORF_count" in row and "microORF" in row
    assert row["cds_extends_beyond_window"] is True
    assert "candidate for extension" in row["notes"]
    assert schema.validate_row(row) == []

def test_sensitivity_bin_helper():
    assert schema.sensitivity_bin_from_delta("") == ""
    assert schema.sensitivity_bin_from_delta(30) == "low"
    assert schema.sensitivity_bin_from_delta(45) == "mid"
    assert schema.sensitivity_bin_from_delta(70) == "high"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(fns)} passed")
    sys.exit(0 if passed == len(fns) else 1)
