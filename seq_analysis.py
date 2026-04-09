"""
Sequence Analysis Engine for PlasmidNet.
Scans a DNA sequence for:
  - Restriction site clusters (cloning artifact signatures)
  - Cryptic promoters (-10 Pribnow box, -35 box)
  - Ribosome Binding Sites (Shine-Dalgarno)
  - Codon usage bias (per gene vs natural background)
  - Natural vs Engineered scoring (k-mer frequency deviation)
"""

import re
import math
from collections import Counter

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
        "engineering_score": 0,  # calculated below
    }

    # Composite engineering score (0-100, higher = more likely engineered)
    results["engineering_score"] = _calculate_engineering_score(results)

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


def _calculate_engineering_score(results):
    """Composite score: higher = more likely engineered."""
    score = 0

    # Restriction site density (natural: ~1-3/kb, engineered: 5+/kb)
    density = results["restriction_density"].get("density_per_kb", 0)
    if density > 8:
        score += 25
    elif density > 5:
        score += 15
    elif density > 3:
        score += 5

    # Restriction site hotspots (clusters = cloning junctions)
    hotspots = len(results["restriction_density"].get("hotspots", []))
    score += min(hotspots * 5, 20)

    # Vector backbone signatures
    vectors = len(results.get("vector_signatures", []))
    score += min(vectors * 15, 30)

    # K-mer naturalness
    kmer = results.get("kmer_score", {})
    score += kmer.get("score", 0) // 3

    # Codon optimization regions
    codon_regions = len(results.get("codon_bias", {}).get("regions", []))
    score += min(codon_regions * 5, 15)

    return min(score, 100)


def _reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(c, c) for c in reversed(seq.upper()))
