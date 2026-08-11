"""Predict regulatory direction (ON / OFF) from the downstream gene annotation.

Rubric (DESIGN.md §2.2):
  * Biosynthesis / salvage genes  -> OFF  (ligand present => feedback repression)
  * Transport / efflux / detox    -> ON   (ligand present => import/expel to rebalance)
  * No product / hypothetical     -> unknown (never guessed)

The class-specific enzyme keywords live here (and only here) so a new metabolic class
is added by extending a keyword list, not by touching the schema. Evidence — the
matched keyword and the gene product string — is always recorded, per task spec.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Generic, class-agnostic keyword rubric ----------------------------------------
OFF_KEYWORDS = [
    "synthase", "synthetase", "biosynthesis", "biosynthetic", "synthesis",
    "reductase", "dehydratase", "cyclase", "kinase", "phosphatase",
    "decarboxylase", "methyltransferase", "aminotransferase", "isomerase",
    "salvage", "oxidoreductase", "oxidase", "dehydrogenase", "pyrophosphorylase",
    "phosphorylase", "deformylase", "hydrolase", "hydroxymethyltransferase",
]
ON_KEYWORDS = [
    "transporter", "permease", "efflux", "importer", "exporter", "import",
    "export", "abc transporter", "abc transport", "transport", "mfs",
    "symporter", "antiporter", "channel", "uptake", "detox", "resistance",
    "extrusion", "porter", "solute binding", "solute-binding",
    "substrate-binding", "substrate binding",
]

# Class-specific overrides / boosters. Keys are Rfam accessions.
# These make specific gene families decisive even when the generic rubric is weak.
CLASS_KEYWORDS = {
    "RF00059": {  # TPP (thiamine pyrophosphate)
        # Biosynthesis + salvage (feedback repression -> OFF).
        # thiO = glycine/amino-acid oxidase; tenA/tenI/thi-4 = thiamine salvage
        # (thiaminase II / aminopyrimidine aminohydrolase); pyrimidine deformylase
        # = HMP-P biosynthesis. Short symbols matched as substrings (see _match).
        "OFF": ["thic", "thid", "thie", "thim", "thig", "thif", "thih", "thio",
                "thiamine biosynthesis", "thiamin biosynthesis",
                "hydroxymethylpyrimidine", "hydroxyethylthiazole", "thiazole",
                "phosphomethylpyrimidine", "methylpyrimidine", "pyrimidine",
                "tena", "teni", "thi-4", "thiaminase", "thiamin-phosphate",
                "thiamine phosphate", "thiamin phosphate"],
        # Thiamine import (ON). tbpA/thiB = periplasmic thiamine-binding protein;
        # ykoEDCF / thiXYZ / pnuT = thiamine ECF/ABC transporters; TonB receptor.
        "ON": ["thit", "thixyz", "thix", "thiy", "thiz", "thib", "tbpa",
               "pnuc", "pnut", "ykoe", "ykod", "ykoc", "ykof",
               "thiamine transporter", "thiamine abc", "thiamine uptake",
               "thiamine-transport", "thiamine transport", "thiamine-binding",
               "thiamine binding", "thiamin-binding", "thiamin binding",
               "tonb-dependent", "siderophore receptor"],
    },
    "RF00162": {  # SAM-I (S box) — methionine/SAM/cysteine biosynthesis & sulfur metab.
        # SAM is the product of Met/SAM biosynthesis -> feedback repression (OFF).
        "OFF": ["meta", "metb", "metc", "mete", "metf", "meth", "metk",
                "metx", "metz", "methionine synthase", "methionine biosynthesis",
                "s-adenosylmethionine synthetase", "sam synthetase", "sam synthase",
                "homoserine", "cystathionine", "cysteine synthase", "cysteine biosynthesis",
                "cysh", "cysk", "cyse", "cysteine desulfurase", "sulfate adenylyltransferase",
                "sulfur", "sulfonate", "sulfite reductase", "trans-sulfuration",
                "o-acetylhomoserine", "aspartate-semialdehyde",
                "adenosylhomocysteinase", "sahh", "spermidine synthase",
                # methionine (SAM) salvage & S-adenosyl-Met production
                "adenosyltransferase", "methylthio", "methionine adenosyltransferase",
                "mtna", "mtnk", "mtnw", "mtnb", "homocysteine methyltransferase"],
        # Methionine / SAM import (ON) — less common for SAM-I.
        "ON": ["metn", "metq", "methionine transporter", "methionine abc",
               "methionine uptake", "methionine import", "metnpq", "d-methionine"],
    },
}


@dataclass
class DirectionCall:
    direction: str = "unknown"       # ON | OFF | unknown
    evidence: str = ""               # matched keyword + product string
    basis: str = "none"             # annotation | architecture | both | none


def _is_gene_symbol(kw: str) -> bool:
    """Short alphanumeric token (e.g. 'meth', 'metk', 'thic', 'tbpa') — must match as a
    whole word so 'meth' (metH) does NOT match inside 'methionine', 'meta' inside
    'metabolism', etc. Multi-word / hyphenated / long descriptive keywords use substring."""
    return len(kw) <= 5 and kw.isalnum()


def _match(text: str, keywords) -> str | None:
    for kw in keywords:
        if _is_gene_symbol(kw):
            # allow one optional trailing digit (metK1, thiC2) but keep the word boundary
            if re.search(r"\b" + re.escape(kw) + r"\d?\b", text):
                return kw
        elif kw in text:
            return kw
    return None


def predict_direction(
    gene_product: str,
    gene_symbol: str = "",
    rfam_family: str = "",
) -> DirectionCall:
    product = (gene_product or "").strip()
    text = f"{product} {gene_symbol}".lower()

    if not product or "hypothetical" in text or product.lower() in {"unknown", "tbd"}:
        return DirectionCall("unknown", f'no informative product ("{product}")', "none")

    # 1) Class-specific keywords take precedence (most decisive).
    cls = CLASS_KEYWORDS.get(rfam_family, {})
    hit = _match(text, cls.get("OFF", []))
    if hit:
        return DirectionCall("OFF", f'{rfam_family} biosynthesis:"{hit}" in "{product}"', "annotation")
    hit = _match(text, cls.get("ON", []))
    if hit:
        return DirectionCall("ON", f'{rfam_family} transport:"{hit}" in "{product}"', "annotation")

    # 2) Generic rubric. Transport wins ties (transport annotations are more specific
    #    than the broad "kinase"/"synthase" verbs that also appear in transport genes).
    on_hit = _match(text, ON_KEYWORDS)
    off_hit = _match(text, OFF_KEYWORDS)
    if on_hit and not off_hit:
        return DirectionCall("ON", f'transport:"{on_hit}" in "{product}"', "annotation")
    if off_hit and not on_hit:
        return DirectionCall("OFF", f'biosynthesis:"{off_hit}" in "{product}"', "annotation")
    if on_hit and off_hit:
        # Ambiguous annotation: prefer transport but say so.
        return DirectionCall(
            "ON",
            f'ambiguous (transport:"{on_hit}" & biosynth:"{off_hit}") in "{product}"; chose ON',
            "annotation",
        )
    return DirectionCall("unknown", f'no rubric keyword in "{product}"', "none")
