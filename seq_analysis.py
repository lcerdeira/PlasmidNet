"""
Sequence Analysis Engine for PlasmidNet.
Scans a DNA sequence for:
  - Restriction site clusters (cloning artifact signatures)
  - Cryptic promoters (-10 Pribnow box, -35 box)
  - Ribosome Binding Sites (Shine-Dalgarno)
  - Codon usage bias (per gene vs natural background)
  - IS elements / transposons (natural mobile element markers)
  - Natural vs Engineered scoring (percentile-based, PLSDB calibrated)
"""

import json
import os
import re
import math
from collections import Counter

# Load PLSDB baseline stats if available
_BASELINE_PATH = os.path.join(os.path.dirname(__file__), "data", "baseline_stats.json")
BASELINE = {}
if os.path.exists(_BASELINE_PATH):
    with open(_BASELINE_PATH) as _f:
        BASELINE = json.load(_f)

# ── Restriction enzyme recognition sequences ───────────────────────
# Common cloning enzymes (Type II, 6+ bp)
RESTRICTION_SITES = {
    # 6-cutters (common cloning)
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
    "SalI": "GTCGAC", "XhoI": "CTCGAG", "NdeI": "CATATG",
    "NcoI": "CCATGG", "XbaI": "TCTAGA", "SpeI": "ACTAGT",
    "PstI": "CTGCAG", "SacI": "GAGCTC", "KpnI": "GGTACC",
    "NotI": "GCGGCCGC", "SfiI": "GGCCNNNNNGGCC",
    "EcoRV": "GATATC", "ClaI": "ATCGAT", "BglII": "AGATCT",
    "NheI": "GCTAGC", "AvrII": "CCTAGG", "MluI": "ACGCGT",
    # Golden Gate / Type IIS
    "BsaI": "GGTCTC", "BbsI": "GAAGAC", "BsmBI": "CGTCTC",
    "SapI": "GCTCTTC", "AarI": "CACCTGC",
    # Gibson Assembly overhangs often leave unique junctions
}

# ── Promoter motifs ────────────────────────────────────────────────
# E. coli sigma-70 consensus
PRIBNOW_BOX = re.compile(r"TA[AT]AAT", re.IGNORECASE)       # -10 region
MINUS35_BOX = re.compile(r"TTGAC[AT]", re.IGNORECASE)        # -35 region
# Spacer between -35 and -10 is typically 15-21 bp
PROMOTER_PATTERN = re.compile(
    r"(TTGAC[AT]).{14,20}(TA[AT]AAT)", re.IGNORECASE
)

# ── Shine-Dalgarno (RBS) ──────────────────────────────────────────
# Consensus: AGGAGG (or partial matches)
SD_SEQUENCES = [
    re.compile(r"AAGGAG", re.IGNORECASE),  # strong
    re.compile(r"AGGAGG", re.IGNORECASE),  # strong
    re.compile(r"AGGAG", re.IGNORECASE),   # moderate
    re.compile(r"GAGG", re.IGNORECASE),    # weak
]

# ── Common vector backbone signatures ──────────────────────────────
VECTOR_SIGNATURES = {
    "pUC_ori": "ATGAGTATTCAACATTTCCGTGTCGCCC",
    "pBR322_ori": "ATTATCATGACATTAACC",
    "ColE1_ori": "ATGAGTATTCAACATTTC",
    "f1_ori": "ACGCTCAGTGGAACGAAAACTCACG",
    "lacZ_alpha": "ATGACCATGATTACGCCAAGC",
    "T7_promoter": "TAATACGACTCACTATAGGG",
    "T3_promoter": "AATTAACCCTCACTAAAGGG",
    "SP6_promoter": "ATTTAGGTGACACTATAG",
    "CMV_promoter": "GACATTGATTATTGACTAG",
    "bla_promoter": "CTCACTCAAAGGCGGTAATAC",
}

# ── IS elements / transposon signatures ────────────────────────────
# Terminal inverted repeats (TIR) of common IS families
IS_SIGNATURES = {
    "IS1": "GGTGATGCTGCCAAC",
    "IS26": "GCACTGTTGCAAATAGTCGGTGGTG",
    "IS903": "TGGATTTATCAGACGATG",
    "IS10": "CTGATGAATCCCCTAATG",
    "IS3": "TGATCAAACTCAAG",
    "IS5": "GAAACAGGCCAGCGG",
    "ISEcp1": "TCTAATCTTGCCTGC",   # mobilises bla genes
    "Tn3_IR": "GGGGAGTGATTTGTTATCATG",
    "Tn21_IR": "GGGGTGTGCTCAAGTAG",
    "Tn1721_IR": "TGATTTTTCATGATCTAG",
}

# Transposase gene fragments (short conserved motifs)
TRANSPOSASE_MOTIFS = [
    re.compile(r"ATGAG[CT]AC[ACGT]AA[CT]GA[CT]GA", re.IGNORECASE),  # DDE domain
    re.compile(r"CC[ACGT]TT[ACGT]GG[ACGT]AA[ACGT]CC", re.IGNORECASE),  # IS helix-turn-helix
]

# Engineering scars are detected by _find_engineering_scars(), not regex
# because they require context (paired sites, specific spacing)


def analyze_sequence(seq, name="query"):
    """
    Run full analysis on a DNA sequence.
    Returns dict with all findings.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    # Remove FASTA header if present
    if seq.startswith(">"):
        lines = seq.split("\n")
        name = lines[0][1:].strip() if lines[0][1:].strip() else name
        seq = "".join(l for l in lines[1:] if not l.startswith(">"))
    seq = re.sub(r"[^ACGTN]", "", seq)

    length = len(seq)
    if length < 100:
        return {"error": "Sequence too short (minimum 100 bp)"}

    results = {
        "name": name,
        "length": length,
        "gc_content": _gc_content(seq),
        "restriction_sites": _find_restriction_sites(seq),
        "restriction_density": _restriction_site_density(seq),
        "promoters": _find_promoters(seq),
        "rbs_sites": _find_rbs(seq),
        "vector_signatures": _find_vector_signatures(seq),
        "codon_bias": _codon_usage_analysis(seq),
        "kmer_score": _kmer_naturalness_score(seq),
        "is_elements": _find_is_elements(seq),
        "engineering_scars": _find_engineering_scars(seq),
        "engineering_score": 0,  # calculated below
        "classification": "",
    }

    # Improved composite score using PLSDB baseline
    results["engineering_score"] = _calculate_engineering_score_v2(results)
    score = results["engineering_score"]
    if score < 20:
        results["classification"] = "Natural"
    elif score < 40:
        results["classification"] = "Natural (resistance platform)"
    elif score < 60:
        results["classification"] = "Ambiguous"
    elif score < 80:
        results["classification"] = "Likely engineered"
    else:
        results["classification"] = "Engineered"

    return results


def _gc_content(seq):
    gc = sum(1 for c in seq if c in "GC")
    return round(gc / len(seq) * 100, 1) if seq else 0


def _find_restriction_sites(seq):
    """Find all restriction enzyme recognition sites and their positions."""
    sites = []
    for enzyme, recognition in RESTRICTION_SITES.items():
        # Handle degenerate bases (N = any)
        pattern = recognition.replace("N", "[ACGT]")
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            sites.append({
                "enzyme": enzyme,
                "position": m.start(),
                "sequence": m.group(),
                "length": len(recognition.replace("N", "")),
            })
        # Also check reverse complement
        rc = _reverse_complement(recognition)
        if rc != recognition:
            rc_pattern = rc.replace("N", "[ACGT]")
            for m in re.finditer(rc_pattern, seq, re.IGNORECASE):
                sites.append({
                    "enzyme": f"{enzyme}(rc)",
                    "position": m.start(),
                    "sequence": m.group(),
                    "length": len(recognition.replace("N", "")),
                })

    return sorted(sites, key=lambda x: x["position"])


def _restriction_site_density(seq):
    """
    Calculate restriction site density in sliding windows.
    Clusters of restriction sites suggest cloning junctions.
    """
    sites = _find_restriction_sites(seq)
    if not sites:
        return {"hotspots": [], "total_sites": 0}

    window = 500
    step = 100
    densities = []
    for start in range(0, len(seq) - window, step):
        end = start + window
        count = sum(1 for s in sites if start <= s["position"] < end)
        if count >= 3:  # hotspot threshold
            densities.append({
                "position": start + window // 2,
                "count": count,
                "start": start, "end": end,
            })

    return {
        "hotspots": densities,
        "total_sites": len(sites),
        "density_per_kb": round(len(sites) / (len(seq) / 1000), 2),
    }


def _find_promoters(seq):
    """Find putative sigma-70 promoter sequences (-35 + spacer + -10)."""
    promoters = []
    for m in PROMOTER_PATTERN.finditer(seq):
        minus35 = m.group(1)
        minus10 = m.group(2)
        spacer_len = m.end(1) - m.start(1)
        # Score based on consensus match
        score = 0
        if minus35.upper() == "TTGACA":
            score += 50
        else:
            score += 30
        if minus10.upper() == "TATAAT":
            score += 50
        else:
            score += 30
        # Optimal spacer is 17 bp
        spacer = m.start(2) - m.end(1)
        if 16 <= spacer <= 18:
            score += 20
        elif 15 <= spacer <= 19:
            score += 10

        promoters.append({
            "position": m.start(),
            "end": m.end(),
            "minus35": minus35,
            "minus10": minus10,
            "spacer": spacer,
            "score": min(score, 100),
            "strand": "+",
        })

    # Also check reverse complement
    rc_seq = _reverse_complement(seq)
    for m in PROMOTER_PATTERN.finditer(rc_seq):
        pos = len(seq) - m.end()
        spacer = m.start(2) - m.end(1)
        score = 60  # simplified for reverse
        if 16 <= spacer <= 18:
            score += 20
        promoters.append({
            "position": pos,
            "end": pos + (m.end() - m.start()),
            "minus35": m.group(1),
            "minus10": m.group(2),
            "spacer": spacer,
            "score": min(score, 100),
            "strand": "-",
        })

    return sorted(promoters, key=lambda x: -x["score"])


def _find_rbs(seq):
    """Find Shine-Dalgarno (ribosome binding site) sequences."""
    rbs_hits = []
    for i, pattern in enumerate(SD_SEQUENCES):
        strength = ["strong", "strong", "moderate", "weak"][i]
        for m in pattern.finditer(seq):
            # Check if an ATG start codon follows within 5-13 bp
            downstream = seq[m.end():m.end() + 15]
            atg_pos = downstream.find("ATG")
            has_atg = atg_pos >= 4 and atg_pos <= 12
            rbs_hits.append({
                "position": m.start(),
                "sequence": m.group(),
                "strength": strength,
                "has_downstream_atg": has_atg,
                "atg_distance": atg_pos if has_atg else None,
            })

    # Deduplicate overlapping hits (keep strongest)
    seen_positions = set()
    unique = []
    for rbs in sorted(rbs_hits, key=lambda x: ["strong", "moderate", "weak"].index(x["strength"])):
        if not any(abs(rbs["position"] - p) < 6 for p in seen_positions):
            unique.append(rbs)
            seen_positions.add(rbs["position"])

    return unique


def _find_vector_signatures(seq):
    """Check for common vector backbone sequences."""
    found = []
    for name, signature in VECTOR_SIGNATURES.items():
        pos = seq.find(signature.upper())
        if pos >= 0:
            found.append({"name": name, "position": pos, "length": len(signature)})
        # Check reverse complement too
        rc = _reverse_complement(signature)
        pos_rc = seq.find(rc.upper())
        if pos_rc >= 0 and pos_rc != pos:
            found.append({"name": f"{name}(rc)", "position": pos_rc, "length": len(signature)})
    return found


def _codon_usage_analysis(seq):
    """
    Analyze codon usage bias across the sequence.
    Compare to E. coli expected frequencies.
    Returns regions with unusual codon bias (potential optimization).
    """
    # E. coli K-12 codon usage (simplified, top codons per amino acid)
    ECOLI_FREQ = {
        "GCG": 0.36, "GCC": 0.27, "GCA": 0.21, "GCT": 0.16,  # Ala
        "CGT": 0.38, "CGC": 0.36, "AGA": 0.04, "AGG": 0.02,  # Arg
        "GAT": 0.63, "GAC": 0.37,  # Asp
        "AAT": 0.45, "AAC": 0.55,  # Asn
        "TGT": 0.44, "TGC": 0.56,  # Cys
        "GAA": 0.69, "GAG": 0.31,  # Glu
        "CAA": 0.35, "CAG": 0.65,  # Gln
        "GGT": 0.34, "GGC": 0.40,  # Gly
        "ATT": 0.51, "ATC": 0.42, "ATA": 0.07,  # Ile
        "CTG": 0.50, "TTA": 0.13, "TTG": 0.13,  # Leu
        "AAA": 0.76, "AAG": 0.24,  # Lys
        "TTT": 0.57, "TTC": 0.43,  # Phe
        "CCG": 0.53, "CCA": 0.19,  # Pro
        "TCT": 0.15, "AGC": 0.28, "TCC": 0.15,  # Ser
        "ACC": 0.44, "ACG": 0.27,  # Thr
        "TAT": 0.57, "TAC": 0.43,  # Tyr
        "GTT": 0.26, "GTC": 0.22, "GTG": 0.37,  # Val
    }

    if len(seq) < 300:
        return {"cai": None, "regions": []}

    # Scan in windows of 300 bp (100 codons)
    window = 300
    step = 150
    regions = []
    for start in range(0, len(seq) - window, step):
        chunk = seq[start:start + window]
        codons = [chunk[i:i+3] for i in range(0, len(chunk) - 2, 3)]
        codons = [c for c in codons if len(c) == 3 and all(b in "ACGT" for b in c)]

        if not codons:
            continue

        # Calculate Codon Adaptation Index (simplified)
        scores = []
        for c in codons:
            if c in ECOLI_FREQ:
                scores.append(ECOLI_FREQ[c])

        if scores:
            cai = math.exp(sum(math.log(max(s, 0.01)) for s in scores) / len(scores))
            if cai > 0.7:  # High CAI = potentially codon-optimized
                regions.append({
                    "start": start, "end": start + window,
                    "cai": round(cai, 3),
                    "note": "High CAI - possible codon optimization",
                })
            elif cai < 0.2:
                regions.append({
                    "start": start, "end": start + window,
                    "cai": round(cai, 3),
                    "note": "Low CAI - unusual codon usage",
                })

    return {"regions": regions}


def _kmer_naturalness_score(seq, k=4):
    """
    Score sequence naturalness based on k-mer frequency deviation
    from expected bacterial plasmid distribution.
    Score 0-100: 0 = very natural, 100 = very engineered.
    """
    if len(seq) < 500:
        return {"score": 50, "note": "Sequence too short for reliable scoring"}

    # Count k-mers
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(c in "ACGT" for c in kmer):
            kmers[kmer] += 1

    total = sum(kmers.values())
    if total == 0:
        return {"score": 50, "note": "No valid k-mers"}

    # Expected frequencies for natural bacterial DNA (uniform-ish with GC bias)
    gc = sum(1 for c in seq if c in "GC") / len(seq)

    # Calculate entropy - natural DNA has higher entropy than engineered
    freqs = [v / total for v in kmers.values()]
    entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
    max_entropy = math.log2(4 ** k)  # maximum for k-mer
    entropy_ratio = entropy / max_entropy

    # Check for palindrome enrichment (restriction sites are palindromic)
    palindromes = 0
    for kmer, count in kmers.items():
        rc = _reverse_complement(kmer)
        if kmer == rc:
            palindromes += count

    palindrome_fraction = palindromes / total if total > 0 else 0

    # Score: low entropy + high palindromes = more likely engineered
    score = 0
    if entropy_ratio < 0.85:
        score += 30
    elif entropy_ratio < 0.90:
        score += 15
    if palindrome_fraction > 0.15:
        score += 25
    elif palindrome_fraction > 0.10:
        score += 10

    # Check for long perfect repeats (common in engineered constructs)
    repeat_score = _check_perfect_repeats(seq)
    score += repeat_score

    return {
        "score": min(score, 100),
        "entropy_ratio": round(entropy_ratio, 3),
        "palindrome_fraction": round(palindrome_fraction, 3),
        "note": "Low" if score < 25 else "Moderate" if score < 50 else "High",
    }


def _check_perfect_repeats(seq, min_len=20):
    """Check for long perfect repeats (common in synthetic constructs)."""
    score = 0
    for length in [50, 30, 20]:
        for i in range(0, len(seq) - length * 2, length):
            segment = seq[i:i + length]
            rest = seq[i + length:]
            if segment in rest:
                if length >= 50:
                    score += 20
                elif length >= 30:
                    score += 10
                else:
                    score += 5
                break
        if score > 0:
            break
    return min(score, 25)


def _find_is_elements(seq):
    """Detect IS element and transposon signatures (natural mobile elements)."""
    found = []
    for name, signature in IS_SIGNATURES.items():
        for m in re.finditer(re.escape(signature), seq, re.IGNORECASE):
            found.append({"name": name, "position": m.start(), "type": "IS/Tn"})
        rc = _reverse_complement(signature)
        for m in re.finditer(re.escape(rc), seq, re.IGNORECASE):
            found.append({"name": f"{name}(rc)", "position": m.start(), "type": "IS/Tn"})

    # Check for transposase motifs
    for pattern in TRANSPOSASE_MOTIFS:
        for m in pattern.finditer(seq):
            found.append({"name": "transposase_motif", "position": m.start(),
                          "type": "transposase"})

    return found


def _find_engineering_scars(seq):
    """
    Detect engineered-specific junction scars.
    Only flags patterns that require PAIRED sites or specific contexts
    that are unlikely in natural DNA.
    """
    found = []

    # Golden Gate: look for PAIRED BsaI sites in convergent orientation
    # (GGTCTC...GAGACC) within 50-5000 bp — this is the assembly signature
    bsai_fwd = [m.start() for m in re.finditer("GGTCTC", seq, re.IGNORECASE)]
    bsai_rev = [m.start() for m in re.finditer("GAGACC", seq, re.IGNORECASE)]
    for f in bsai_fwd:
        for r in bsai_rev:
            dist = abs(r - f)
            if 50 < dist < 5000:
                found.append({"name": "Golden_Gate_pair",
                              "position": min(f, r),
                              "sequence": f"BsaI pair at {f:,} & {r:,} ({dist} bp apart)"})

    # Multiple Cloning Site (MCS): 4+ different 6-cutter sites within 200 bp
    sites_by_pos = {}
    for enzyme, recognition in RESTRICTION_SITES.items():
        if len(recognition.replace("N", "")) >= 6:
            pattern = recognition.replace("N", "[ACGT]")
            for m in re.finditer(pattern, seq, re.IGNORECASE):
                sites_by_pos.setdefault(m.start() // 200 * 200, set()).add(enzyme)
    for window_start, enzymes in sites_by_pos.items():
        if len(enzymes) >= 4:
            found.append({"name": "Multiple_Cloning_Site",
                          "position": window_start,
                          "sequence": f"{len(enzymes)} unique 6-cutters in 200bp: {', '.join(sorted(enzymes)[:5])}"})

    return found


def _re_hotspots_near_is(results):
    """
    Check if RE hotspots coincide with IS element positions.
    If they do, the hotspots are likely natural transposon boundaries,
    not cloning junctions.
    """
    is_positions = set()
    for is_el in results.get("is_elements", []):
        # Extend IS position to a window
        pos = is_el["position"]
        for p in range(max(0, pos - 500), pos + 500):
            is_positions.add(p)

    hotspots = results.get("restriction_density", {}).get("hotspots", [])
    near_is = 0
    for hs in hotspots:
        if hs["position"] in is_positions:
            near_is += 1

    return near_is, len(hotspots)


def _calculate_engineering_score_v2(results):
    """
    Improved engineering score using PLSDB baseline calibration.
    Key improvements:
      - RE hotspots near IS elements are discounted (natural transposon boundaries)
      - Vector signatures are weighted by specificity (ColE1 is also natural)
      - Codon bias is species-aware
      - Engineering scars (Golden Gate/Gibson) are strong indicators
    """
    score = 0
    reasons = []

    # 1. Engineering-specific scars (strongest signal)
    eng_scars = len(results.get("engineering_scars", []))
    if eng_scars > 0:
        score += min(eng_scars * 20, 40)
        reasons.append(f"{eng_scars} engineering junction scars")

    # 2. Vector signatures — but discount shared natural/synthetic ones
    vectors = results.get("vector_signatures", [])
    synthetic_only = {"T7_promoter", "T3_promoter", "SP6_promoter",
                      "CMV_promoter", "lacZ_alpha", "f1_ori"}
    shared_natural = {"pBR322_ori", "ColE1_ori", "pUC_ori", "bla_promoter"}
    for v in vectors:
        base_name = v["name"].replace("(rc)", "")
        if base_name in synthetic_only:
            score += 15  # strong engineered signal
            reasons.append(f"{base_name} (synthetic-specific)")
        elif base_name in shared_natural:
            score += 3   # weak signal — also found in natural IncQ/ColE1 plasmids
            reasons.append(f"{base_name} (shared natural/synthetic)")

    # 3. RE density — compare to PLSDB baseline
    density = results["restriction_density"].get("density_per_kb", 0)
    nat_p95 = BASELINE.get("re_density_natural_p95", 6.0)
    eng_thresh = BASELINE.get("re_density_engineered_threshold", 8.0)
    if density > eng_thresh:
        score += 15
    elif density > nat_p95:
        score += 5

    # 4. RE hotspots — discount those near IS elements
    near_is, total_hs = _re_hotspots_near_is(results)
    natural_hs = near_is  # hotspots explained by transposons
    unexplained_hs = total_hs - natural_hs
    if unexplained_hs > 5:
        score += 15
    elif unexplained_hs > 2:
        score += 5
    # Bonus for IS elements (natural indicators, REDUCE score)
    is_count = len(results.get("is_elements", []))
    if is_count >= 3:
        score -= 15
        reasons.append(f"{is_count} IS elements (natural mobile element markers)")
    elif is_count >= 1:
        score -= 5

    # 5. K-mer naturalness
    kmer = results.get("kmer_score", {})
    kmer_s = kmer.get("score", 0)
    score += kmer_s // 5  # reduced weight vs v1

    # 6. Codon optimization (only if strong signal)
    codon_regions = results.get("codon_bias", {}).get("regions", [])
    high_cai = [r for r in codon_regions if r.get("cai", 0) > 0.8]
    if len(high_cai) >= 3:
        score += 10

    # Floor at 0
    score = max(0, min(score, 100))

    results["score_reasons"] = reasons
    return score


def _reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(c, c) for c in reversed(seq.upper()))
