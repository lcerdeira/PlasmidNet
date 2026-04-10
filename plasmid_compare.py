"""
Plasmid comparison: align and visualize up to 10 plasmids.
Uses pyCirclize for publication-quality circular comparison plots.
"""

import base64
import io
import os
import re
import subprocess
import tempfile


def compare_plasmids(accessions, sequences, names=None):
    """
    Compare multiple plasmid sequences.
    Returns a base64-encoded PNG of the comparison plot + alignment data.

    Args:
        accessions: list of accession strings
        sequences: list of DNA sequence strings
        names: optional list of display names
    """
    if not names:
        names = accessions

    n = len(sequences)
    if n < 2 or n > 10:
        return None, "Provide 2-10 plasmids for comparison."

    # Run pairwise BLASTn alignments
    alignments = _run_pairwise_blast(sequences, names)

    # Generate comparative circular plot
    img_b64 = _make_comparison_plot(sequences, names, alignments)

    return img_b64, alignments


def _run_pairwise_blast(sequences, names):
    """Run BLASTn between all pairs. Returns alignment results."""
    alignments = []

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write all sequences as FASTA
        fasta_files = []
        for i, (seq, name) in enumerate(zip(sequences, names)):
            path = os.path.join(tmpdir, f"seq_{i}.fasta")
            with open(path, "w") as f:
                f.write(f">{name}\n{seq}\n")
            fasta_files.append(path)

        # Make BLAST database from the first sequence and align others
        ref = fasta_files[0]
        db_path = os.path.join(tmpdir, "refdb")

        try:
            subprocess.run(
                ["makeblastdb", "-in", ref, "-dbtype", "nucl", "-out", db_path],
                capture_output=True, timeout=30,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return []

        for i in range(1, len(fasta_files)):
            try:
                result = subprocess.run(
                    ["blastn", "-query", fasta_files[i], "-db", db_path,
                     "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue",
                     "-max_hsps", "50", "-evalue", "1e-10"],
                    capture_output=True, text=True, timeout=60,
                )
                for line in result.stdout.strip().split("\n"):
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 9:
                        alignments.append({
                            "query": names[i],
                            "subject": names[0],
                            "identity": float(parts[2]),
                            "length": int(parts[3]),
                            "q_start": int(parts[4]),
                            "q_end": int(parts[5]),
                            "s_start": int(parts[6]),
                            "s_end": int(parts[7]),
                            "evalue": float(parts[8]),
                        })
            except (subprocess.TimeoutExpired, Exception):
                continue

    return alignments


def _make_comparison_plot(sequences, names, alignments):
    """Create a pyCirclize multi-ring comparison plot."""
    import matplotlib
    matplotlib.use("Agg")
    from pycirclize import Circos
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # Build sector sizes
    sectors = {name: len(seq) for name, seq in zip(names, sequences)}

    # Create Circos with one sector per plasmid
    circos = Circos(sectors, space=8)

    # Color palette for tracks
    colors = ["#4a90d9", "#5cb85c", "#e74c3c", "#f0ad4e", "#9b59b6",
              "#1abc9c", "#e67e22", "#3498db", "#e91e63", "#00bcd4"]

    for i, sector in enumerate(circos.sectors):
        color = colors[i % len(colors)]
        seq = sequences[i]

        # Outer track: sequence backbone
        track = sector.add_track((85, 95))
        track.axis(fc=color, alpha=0.3, ec=color, lw=0.5)

        # GC content inner track
        gc_track = sector.add_track((70, 83))
        gc_track.axis(fc="#fafafa", ec="#e2ddd5", lw=0.3)

        # Compute GC in windows
        window = max(len(seq) // 100, 100)
        positions = []
        gc_vals = []
        for start in range(0, len(seq) - window, window):
            chunk = seq[start:start + window]
            gc = sum(1 for c in chunk if c in "GC") / len(chunk) * 100
            positions.append(start + window // 2)
            gc_vals.append(gc)

        if positions:
            gc_mean = sum(gc_vals) / len(gc_vals)
            gc_track.fill_between(
                positions, gc_vals, [gc_mean] * len(positions),
                fc=color, alpha=0.3, ec="none",
            )
            gc_track.line(positions, gc_vals, color=color, lw=0.5)

        # Add tick marks
        tick_interval = max(len(seq) // 5, 1000)
        ticks = list(range(0, len(seq), tick_interval))
        labels = [f"{t // 1000}kb" if t >= 1000 else str(t) for t in ticks]
        if labels:
            labels[0] = "0"
        track.xticks(ticks, labels, label_size=6, label_margin=1)

    # Draw alignment links between sectors
    if alignments and len(circos.sectors) >= 2:
        ref_sector = circos.sectors[0]
        for aln in alignments:
            # Find query sector
            query_sector = None
            for s in circos.sectors:
                if s.name == aln["query"]:
                    query_sector = s
                    break
            if not query_sector:
                continue

            identity = aln["identity"]
            if identity >= 99:
                link_color = "#2ecc71"
                alpha = 0.4
            elif identity >= 95:
                link_color = "#3498db"
                alpha = 0.3
            elif identity >= 90:
                link_color = "#f39c12"
                alpha = 0.25
            else:
                link_color = "#e74c3c"
                alpha = 0.2

            circos.link(
                (ref_sector.name, aln["s_start"], aln["s_end"]),
                (query_sector.name, aln["q_start"], aln["q_end"]),
                color=link_color, alpha=alpha,
            )

    # Center text
    circos.text(f"Comparison\n{len(names)} plasmids", size=8, r=0)

    # Render
    fig = circos.plotfig(dpi=150)

    # Legend
    legend_items = [
        mpatches.Patch(color=colors[i % len(colors)], alpha=0.3,
                       label=f"{names[i]} ({len(sequences[i]):,} bp)")
        for i in range(len(names))
    ]
    if alignments:
        legend_items.extend([
            mpatches.Patch(color="#2ecc71", alpha=0.4, label=">=99% identity"),
            mpatches.Patch(color="#3498db", alpha=0.3, label="95-99%"),
            mpatches.Patch(color="#f39c12", alpha=0.25, label="90-95%"),
            mpatches.Patch(color="#e74c3c", alpha=0.2, label="<90%"),
        ])
    fig.legend(handles=legend_items, loc="lower center",
               ncol=min(4, len(legend_items)), fontsize=6, frameon=False)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white", dpi=150)
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def export_genbank(accession, sequence, features=None, organism="Unknown"):
    """Generate a GenBank format string for a plasmid."""
    length = len(sequence)
    gc = sum(1 for c in sequence if c in "GC") / length * 100

    lines = []
    lines.append(f"LOCUS       {accession:16s} {length:>8d} bp    DNA     circular     UNK")
    lines.append(f"DEFINITION  {accession} plasmid, complete sequence.")
    lines.append(f"ACCESSION   {accession}")
    lines.append(f"SOURCE      {organism}")
    lines.append(f"  ORGANISM  {organism}")
    lines.append("FEATURES             Location/Qualifiers")
    lines.append(f"     source          1..{length}")
    lines.append(f'                     /organism="{organism}"')
    lines.append(f'                     /mol_type="genomic DNA"')

    if features:
        for f in features:
            start = f.get("start", 1)
            end = f.get("end", length)
            strand = f.get("strand", 1)
            gene = f.get("gene", "")
            product = f.get("product", "")
            if strand == -1:
                loc = f"complement({start}..{end})"
            else:
                loc = f"{start}..{end}"
            lines.append(f"     CDS             {loc}")
            if gene:
                lines.append(f'                     /gene="{gene}"')
            if product:
                lines.append(f'                     /product="{product}"')

    lines.append("ORIGIN")
    for i in range(0, length, 60):
        chunk = sequence[i:i + 60].lower()
        parts = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9d} {parts}")
    lines.append("//")

    return "\n".join(lines)


def export_fasta(accession, sequence, description=""):
    """Generate FASTA format string."""
    lines = [f">{accession} {description}"]
    for i in range(0, len(sequence), 70):
        lines.append(sequence[i:i + 70])
    return "\n".join(lines)
