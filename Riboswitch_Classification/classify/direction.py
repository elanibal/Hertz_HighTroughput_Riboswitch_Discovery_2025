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

from dataclasses import dataclass

# Generic, class-agnostic keyword rubric ----------------------------------------
OFF_KEYWORDS = [
    "synthase", "synthetase", "biosynthesis", "biosynthetic", "synthesis",
    "reductase", "dehydratase", "cyclase", "kinase", "phosphatase",
    "decarboxylase", "methyltransferase", "aminotransferase", "isomerase",
    "salvage",
]
ON_KEYWORDS = [
    "transporter", "permease", "efflux", "importer", "exporter", "import",
    "export", "abc transporter", "mfs", "symporter", "antiporter", "channel",
    "uptake", "detox", "resistance", "extrusion", "porter",
]

# Class-specific overrides / boosters. Keys are Rfam accessions.
# These make specific gene families decisive even when the generic rubric is weak.
CLASS_KEYWORDS = {
    "RF00059": {  # TPP (thiamine pyrophosphate)
        "OFF": ["thic", "thid", "thie", "thim", "thig", "thif", "thih",
                "thiamine biosynthesis", "hydroxymethylpyrimidine",
                "thiazole", "thiamin biosynthesis", "phosphomethylpyrimidine"],
        "ON": ["thit", "thixyz", "pnuc", "ykoe", "ykod", "ykoc",
               "thiamine transporter", "thiamine abc", "thiamine uptake",
               "thiamine-transport"],
    },
}


@dataclass
class DirectionCall:
    direction: str = "unknown"       # ON | OFF | unknown
    evidence: str = ""               # matched keyword + product string
    basis: str = "none"             # annotation | architecture | both | none


def _match(text: str, keywords) -> str | None:
    for kw in keywords:
        if kw in text:
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
