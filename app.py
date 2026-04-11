"""
PlasmidNet Dashboard - Interactive explorer for the PLSDB 2025 plasmid database.
Built with Dash + Plotly. Deployable on Heroku / AWS free tier.
"""

import math
import os
import logging

import dash
from dash import dcc, html, dash_table, callback, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from data_loader import fetch_genbank_full, _fetch_genbank_text, _parse_genbank_sequence
import db
import seq_analysis
import plasmid_compare

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load data (all from local SQLite — instant, full database)
# ---------------------------------------------------------------------------
logger.info("Loading data from SQLite ...")
overview = db.overview_stats()
taxonomy = db.top_genera(20)
amr = db.amr_gene_counts(30)
amr_with_class = db.amr_genes_with_class(30)
amr_classes = db.amr_drug_class_counts()
mob_typing = {
    "rep_types": db.inc_group_counts(20),
    "relaxase_types": db.relaxase_counts(),
    "mpf_types": db.mpf_counts(),
    "mobility": db.mobility_distribution(),
}
import json as _json
_corr_cache = os.path.join(os.path.dirname(__file__), "data", "correlations_cache.json")
if os.path.exists(_corr_cache):
    with open(_corr_cache) as _f:
        correlations = _json.load(_f)
    logger.info("Loaded correlations from cache")
else:
    correlations = db.build_correlations()
    with open(_corr_cache, "w") as _f:
        _json.dump(correlations, _f)
    logger.info("Built and cached correlations")
temporal = db.temporal_distribution()
gc_dist = db.gc_distribution()
length_dist = db.length_distribution()
logger.info(f"Loaded {overview['total']:,} plasmids from SQLite")

# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------
app = dash.Dash(
    __name__,
    title="PlasmidNet - Dashboard",
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    suppress_callback_exceptions=True,
)
server = app.server  # for gunicorn

# Register REST API
from api import api_bp
server.register_blueprint(api_bp)

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
COLORS = {
    "bg": "#faf7f2",
    "card": "#ffffff",
    "card_border": "#e2ddd5",
    "accent": "#0d9488",
    "accent2": "#6366f1",
    "accent3": "#10b981",
    "accent4": "#f97316",
    "accent5": "#ec4899",
    "text": "#1e293b",
    "text_muted": "#64748b",
    "danger": "#ef4444",
}

PLOTLY_TEMPLATE = "plotly_white"
CHART_COLORS = px.colors.qualitative.Set2


# ---------------------------------------------------------------------------
# Helper: stat card
# ---------------------------------------------------------------------------
def _parse_multi_fasta(text):
    """Parse FASTA text (possibly multi-record) into list of (name, sequence)."""
    import re
    results = []
    text = text.strip()
    if not text.startswith(">"):
        # Raw sequence
        seq = re.sub(r"[^ACGTNacgtn]", "", text)
        if seq:
            results.append(("uploaded", seq.upper()))
        return results
    current_name = ""
    current_seq = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_name and current_seq:
                results.append((current_name, "".join(current_seq).upper()))
            current_name = line[1:].split()[0] or "unnamed"
            current_seq = []
        else:
            current_seq.append(re.sub(r"[^ACGTNacgtn]", "", line))
    if current_name and current_seq:
        results.append((current_name, "".join(current_seq).upper()))
    return results


def stat_card(title, value, icon, color):
    return html.Div(
        className="stat-card",
        children=[
            html.Div(className="stat-icon", children=icon, style={"color": color}),
            html.Div([
                html.H3(f"{value:,}" if isinstance(value, (int, float)) else value,
                         className="stat-value"),
                html.P(title, className="stat-label"),
            ]),
        ],
    )


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------

def make_topology_chart():
    labels = ["Circular", "Linear"]
    values = [overview.get("topology_circular", 0), overview.get("topology_linear", 0)]
    fig = go.Figure(go.Pie(
        labels=labels, values=values,
        hole=0.55, marker_colors=[COLORS["accent"], COLORS["accent4"]],
        textinfo="label+percent", textfont_size=13,
        hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>Percent: %{percent}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        showlegend=False,
    )
    return fig


def make_source_chart():
    labels = ["RefSeq", "INSDC"]
    values = [overview.get("source_RefSeq", 0), overview.get("source_INSDC", 0)]
    fig = go.Figure(go.Pie(
        labels=labels, values=values,
        hole=0.55, marker_colors=[COLORS["accent2"], COLORS["accent3"]],
        textinfo="label+percent", textfont_size=13,
        hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>Percent: %{percent}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        showlegend=False,
    )
    return fig


def make_taxonomy_chart():
    if not taxonomy:
        return go.Figure()
    df = pd.DataFrame(
        sorted(taxonomy.items(), key=lambda x: x[1], reverse=True),
        columns=["Genus", "Plasmid Count"],
    )
    fig = px.bar(
        df, x="Plasmid Count", y="Genus", orientation="h",
        color="Plasmid Count", color_continuous_scale="Viridis",
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=550,
        yaxis=dict(autorange="reversed"), coloraxis_showscale=False,
        xaxis_title="Number of Plasmids", yaxis_title="",
    )
    return fig


def make_amr_chart():
    """Bar chart of top AMR genes coloured by their real drug class from the DB."""
    if not amr_with_class:
        return go.Figure()
    df = pd.DataFrame(amr_with_class)
    # Normalise class names for display
    class_map = {
        "BETA-LACTAM": "Beta-lactam", "AMINOGLYCOSIDE": "Aminoglycoside",
        "aminoglycoside antibiotic": "Aminoglycoside",
        "SULFONAMIDE": "Sulfonamide", "sulfonamide antibiotic": "Sulfonamide",
        "TETRACYCLINE": "Tetracycline", "tetracycline antibiotic": "Tetracycline",
        "QUATERNARY AMMONIUM": "QAC", "PHENICOL": "Phenicol",
        "macrolide antibiotic": "Macrolide", "TRIMETHOPRIM": "Trimethoprim",
        "MERCURY": "Mercury", "COPPER": "Copper", "COPPER/SILVER": "Copper/Silver",
        "TELLURIUM": "Tellurium", "ARSENIC": "Arsenic",
        "QUINOLONE": "Quinolone", "fluoroquinolone antibiotic": "Quinolone",
        "FOSFOMYCIN": "Fosfomycin", "COLISTIN": "Colistin",
        "GLYCOPEPTIDE": "Glycopeptide",
    }
    df["Drug Class"] = df["drug_class"].map(class_map).fillna(df["drug_class"])
    fig = px.bar(
        df, x="gene", y="count", color="Drug Class",
        color_discrete_sequence=px.colors.qualitative.Set2 + px.colors.qualitative.Pastel,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        xaxis_title="", yaxis_title="Number of Plasmids",
        xaxis_tickangle=-45, legend=dict(
            orientation="h", yanchor="bottom", y=1.02,
            xanchor="center", x=0.5, font_size=9,
        ),
    )
    return fig


def make_temporal_chart():
    if not temporal:
        return go.Figure()
    df = pd.DataFrame(
        sorted(temporal.items(), key=lambda x: x[0]),
        columns=["Year", "Count"],
    )
    fig = px.bar(
        df, x="Year", y="Count",
        color_discrete_sequence=[COLORS["accent"]],
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        xaxis_title="Creation Year", yaxis_title="Number of Plasmids",
    )
    return fig


def make_gc_chart():
    if not gc_dist:
        return go.Figure()
    df = pd.DataFrame(
        list(gc_dist.items()), columns=["GC Range (%)", "Count"],
    )
    fig = px.bar(
        df, x="GC Range (%)", y="Count",
        color_discrete_sequence=[COLORS["accent3"]],
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        xaxis_title="GC Content (%)", yaxis_title="Number of Plasmids",
    )
    return fig


def make_length_chart():
    if not length_dist:
        return go.Figure()
    df = pd.DataFrame(
        list(length_dist.items()), columns=["Size Range", "Count"],
    )
    fig = px.bar(
        df, x="Size Range", y="Count",
        color_discrete_sequence=[COLORS["accent2"]],
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        xaxis_title="Plasmid Length", yaxis_title="Number of Plasmids",
    )
    return fig


def make_kingdom_chart():
    labels = ["Bacteria", "Archaea"]
    values = [overview.get("kingdom_Bacteria", 0), overview.get("kingdom_Archaea", 0)]
    fig = go.Figure(go.Pie(
        labels=labels, values=values,
        hole=0.55, marker_colors=[COLORS["accent3"], COLORS["accent5"]],
        textinfo="label+percent", textfont_size=13,
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300,
        showlegend=False,
    )
    return fig


def make_drug_class_summary():
    """Bar chart of AMR drug class distribution from full database."""
    if not amr_classes:
        return go.Figure()
    # Filter out empty class, heavy metals (shown elsewhere), take top 20
    skip = {"", "COPPER", "COPPER/SILVER", "MERCURY", "SILVER",
            "TELLURIUM", "ARSENIC", "ZINC", "LEAD", "CADMIUM", "CHROMIUM"}
    filtered = {k: v for k, v in amr_classes.items() if k.upper() not in skip}
    df = pd.DataFrame(
        sorted(filtered.items(), key=lambda x: x[1], reverse=True)[:20],
        columns=["Drug Class", "Plasmid Count"],
    )
    fig = px.bar(
        df, x="Plasmid Count", y="Drug Class", orientation="h",
        color="Drug Class", color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=500,
        yaxis=dict(autorange="reversed"), showlegend=False,
        xaxis_title="Plasmids with Resistance", yaxis_title="",
    )
    return fig


# ---------------------------------------------------------------------------
# Circular Plasmid Architecture Map  (interactive Plotly polar)
# ---------------------------------------------------------------------------

FEAT_COLORS = {
    "CDS_fwd": "#4a90d9",
    "CDS_rev": "#5cb85c",
    "AMR": "#e74c3c",
    "MOB": "#f0ad4e",
    "BGC": "#9b59b6",
    "TA": "#e67e22",
    "VIR": "#c0392b",
    "QAC": "#1abc9c",
    "GC_high": "#0d9488",
    "GC_low": "#e74c3c",
    "GC_mean": "#94a3b8",
}


def make_plasmid_map(plasmid_features: dict, genbank_features: list = None,
                     gc_profile: list = None):
    """
    Draw an interactive circular plasmid map using Plotly polar charts.
    Returns a go.Figure with zoom and hover-tooltip support.

    Rings (outside-in):
      Forward CDS | Reverse CDS | Classified features | GC content (Barpolar)
    Genes are colour-coded by classify_gene(): AMR, MOB, TA, VIR, QAC, Metal.
    """
    fig = go.Figure()
    length = plasmid_features.get("length", 1)
    if length == 0:
        length = 1

    def bp_to_deg(bp):
        return (bp / length) * 360

    # ── Ring radii (base, bar-height) ───────────────────────────────
    R_FWD  = (9.0, 1.5)   # forward CDS
    R_REV  = (7.0, 1.5)   # reverse CDS
    R_FEAT = (5.0, 1.5)   # PLSDB features (AMR/MOB/BGC with positions)
    R_GC   = (2.5, 2.0)   # GC content bars

    # Map classify_gene() categories → colour
    CAT_COLOR = {
        "AMR":   FEAT_COLORS["AMR"],
        "MOB":   FEAT_COLORS["MOB"],
        "TA":    FEAT_COLORS["TA"],
        "VIR":   FEAT_COLORS["VIR"],
        "QAC":   FEAT_COLORS["QAC"],
        "Metal": "#8e44ad",
    }

    # ── helper: add a Barpolar ring ─────────────────────────────────
    def _add_ring(feats, base, height, default_color, name, lg,
                  hover_fn, show_legend=True):
        if not feats:
            return
        thetas, widths, hovers, colors_list = [], [], [], []
        for f in feats:
            s, e = f["start"], f["end"]
            thetas.append(bp_to_deg((s + e) / 2))
            widths.append(max(bp_to_deg(e - s), 0.6))
            hovers.append(hover_fn(f))
            colors_list.append(f.get("_color", default_color))
        fig.add_trace(go.Barpolar(
            r=[height] * len(feats), theta=thetas, width=widths, base=base,
            marker_color=colors_list, marker_line_color="white",
            marker_line_width=0.4, opacity=0.88,
            name=name, legendgroup=lg, showlegend=show_legend,
            hovertext=hovers, hoverinfo="text",
        ))

    # ── helper: CDS hover text ──────────────────────────────────────
    def _cds_hover(f, strand_label):
        gene = f.get("gene") or f.get("locus_tag") or ""
        cat = f.get("_cat", "CDS")
        cat_line = f"<b>[{cat}]</b> " if cat != "CDS" else ""
        return (
            f"{cat_line}<b>{gene}</b><br>"
            f"{f.get('product','')}<br>"
            f"{f['start']:,} – {f['end']:,} bp<br>"
            f"Strand: {strand_label}"
        )

    # ── Classify GenBank CDS genes ──────────────────────────────────
    if genbank_features:
        fwd, rev = [], []
        for f in genbank_features:
            gene = f.get("gene") or f.get("locus_tag") or ""
            product = f.get("product", "")
            cat = classify_gene(gene, product)
            color = CAT_COLOR.get(cat, FEAT_COLORS["CDS_fwd"] if f.get("strand", 1) == 1
                                  else FEAT_COLORS["CDS_rev"])
            entry = {**f, "_color": color, "_cat": cat}
            if f.get("strand", 1) == 1:
                fwd.append(entry)
            else:
                rev.append(entry)

        _add_ring(fwd, *R_FWD, FEAT_COLORS["CDS_fwd"], "CDS forward", "cds_fwd",
                  lambda f: _cds_hover(f, "forward (+)"))
        _add_ring(rev, *R_REV, FEAT_COLORS["CDS_rev"], "CDS reverse", "cds_rev",
                  lambda f: _cds_hover(f, "reverse (–)"))

    # ── PLSDB positioned features (AMR / MOB / BGC) ────────────────
    for feats, color, name, lg, hover_fn in [
        (plasmid_features.get("amr", []), FEAT_COLORS["AMR"], "AMR (PLSDB)", "amr_p",
         lambda f: f"<b>AMR: {f.get('gene','')}</b><br>{f.get('product','')}<br>"
                   f"Drug: {f.get('drug_class','') or f.get('agent','')}<br>"
                   f"{f['start']:,} – {f['end']:,} bp"),
        (plasmid_features.get("mob", []), FEAT_COLORS["MOB"], "MOB (PLSDB)", "mob_p",
         lambda f: f"<b>MOB: {f.get('element','')}</b><br>"
                   f"Biomarker: {f.get('biomarker','')}<br>"
                   f"{f['start']:,} – {f['end']:,} bp"),
        (plasmid_features.get("bgc", []), FEAT_COLORS["BGC"], "BGC (PLSDB)", "bgc_p",
         lambda f: f"<b>BGC: {f.get('type','')}</b><br>"
                   f"Completeness: {f.get('completeness','')}<br>"
                   f"{f['start']:,} – {f['end']:,} bp"),
    ]:
        _add_ring(feats, *R_FEAT, color, name, lg, hover_fn)

    # ── GC content ring (Barpolar, coloured above/below mean) ──────
    if gc_profile:
        gc_vals = [p["gc"] for p in gc_profile]
        gc_mean = sum(gc_vals) / len(gc_vals)
        gc_max_dev = max(abs(v - gc_mean) for v in gc_vals) or 0.01

        gc_thetas, gc_widths, gc_heights, gc_colors, gc_hovers = [], [], [], [], []
        step_deg = bp_to_deg(gc_profile[1]["pos"] - gc_profile[0]["pos"]) if len(gc_profile) > 1 else 1
        for p in gc_profile:
            gc_thetas.append(bp_to_deg(p["pos"]))
            gc_widths.append(step_deg)
            dev = (p["gc"] - gc_mean) / gc_max_dev
            gc_heights.append(R_GC[1] / 2 + dev * R_GC[1] / 2)
            gc_colors.append(FEAT_COLORS["GC_high"] if p["gc"] >= gc_mean
                             else FEAT_COLORS["GC_low"])
            gc_hovers.append(f"GC: {p['gc']*100:.1f}%  (mean {gc_mean*100:.1f}%)")

        fig.add_trace(go.Barpolar(
            r=gc_heights, theta=gc_thetas, width=gc_widths,
            base=R_GC[0],
            marker_color=gc_colors, marker_line_width=0, opacity=0.6,
            name="GC content", legendgroup="gc",
            hovertext=gc_hovers, hoverinfo="text",
        ))

    # ── Gene name annotations ───────────────────────────────────────
    if genbank_features:
        for f in genbank_features:
            gene = f.get("gene") or ""
            if not gene:
                continue
            mid_deg = bp_to_deg((f["start"] + f["end"]) / 2)
            strand = f.get("strand", 1)
            r_lbl = R_FWD[0] + R_FWD[1] + 0.5 if strand == 1 else R_REV[0] - 0.5
            angle_rad = math.radians(90 - mid_deg)
            fig.add_annotation(
                x=r_lbl * math.cos(angle_rad),
                y=r_lbl * math.sin(angle_rad),
                text=f"<b>{gene}</b>", showarrow=False,
                font=dict(size=8, color=COLORS["text"]),
            )

    # ── Tick labels & center text ───────────────────────────────────
    tick_interval = _nice_tick_interval(length)
    tick_vals = list(range(0, length, tick_interval))
    tick_degs = [bp_to_deg(t) for t in tick_vals]
    tick_labels = [f"{t//1000} kb" if t >= 1000 else str(t) for t in tick_vals]
    tick_labels[0] = "0"

    acc = plasmid_features.get("accession", "")
    topo = plasmid_features.get("topology", "circular")
    gc_val = plasmid_features.get("gc", 0)
    gc_pct = gc_val * 100 if gc_val <= 1 else gc_val
    organism = plasmid_features.get("organism", "")
    if organism:
        organism = organism.split(" (")[0].replace("_", " ")
    center = f"<b>{acc}</b><br>{_format_size(length)} | {topo}"
    if gc_pct > 0:
        center += f" | GC {gc_pct:.1f}%"
    if organism:
        center += f"<br><i>{organism}</i>"

    fig.add_annotation(
        x=0.5, y=0.5, xref="paper", yref="paper",
        text=center, showarrow=False,
        font=dict(size=11, color=COLORS["text"]), align="center",
    )

    # ── Layout ──────────────────────────────────────────────────────
    outer = R_FWD[0] + R_FWD[1] + 2
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=30, b=30, l=30, r=30), height=720,
        polar=dict(
            radialaxis=dict(visible=False, range=[0, outer]),
            angularaxis=dict(
                direction="clockwise", rotation=90,
                tickvals=tick_degs, ticktext=tick_labels,
                tickfont=dict(size=9, color=COLORS["text_muted"]),
                gridcolor="#e2ddd5", linecolor="#cbd5e1",
            ),
            bgcolor="rgba(0,0,0,0)",
        ),
        legend=dict(
            orientation="h", yanchor="top", y=-0.03,
            xanchor="center", x=0.5, font_size=10,
            bgcolor="rgba(255,255,255,0.85)",
        ),
        hoverlabel=dict(
            bgcolor="#ffffff", bordercolor="#e2ddd5",
            font_size=12, font_color="#1e293b",
        ),
        dragmode="zoom",
    )
    return fig


def _is_amr_gene(gene_name: str, product: str) -> bool:
    """Check if a gene is an AMR gene based on name or product."""
    amr_prefixes = ("bla", "aph", "aac", "ant", "aad", "mec", "van", "erm",
                    "tet", "sul", "dfr", "qnr", "mcr", "fos", "cat", "cfr",
                    "mph", "lin", "str", "flo")
    g = gene_name.lower()
    if any(g.startswith(p) for p in amr_prefixes):
        return True
    p = product.lower()
    return any(kw in p for kw in ("beta-lactamase", "carbapenem",
                                   "aminoglycoside", "phosphotransferase"))


def _is_mob_gene(gene_name: str, product: str) -> bool:
    """Check if a gene is a mobilization / conjugation gene."""
    mob_prefixes = ("mob", "tra", "trb", "trc", "trf", "trh", "tri", "trk",
                    "trn", "tro", "trp", "tru", "trw", "virb", "vird")
    g = gene_name.lower()
    if any(g.startswith(p) for p in mob_prefixes):
        return True
    p = product.lower()
    return any(kw in p for kw in ("mobilization", "conjugal", "conjugative",
                                   "type iv secretion", "mating pair",
                                   "coupling protein", "relaxase"))


def _is_ta_gene(gene_name: str, product: str) -> bool:
    """Check if a gene is a toxin-antitoxin system component."""
    ta_names = {"rele", "relb", "higa", "higb", "mazf", "maze", "ccda", "ccdb",
                "pard", "pare", "parc", "hicb", "hica", "vapb", "vapc",
                "phd", "doc", "yoeb", "yefm", "mqsr", "mqsa", "hipb", "hipa",
                "pema", "pemk", "zeta", "epsilon", "fic", "yafq", "dina",
                "brna", "brnb", "pasa", "pasb", "pasc", "rnla", "rnlb"}
    if gene_name.lower() in ta_names:
        return True
    p = product.lower()
    return any(kw in p for kw in ("toxin-antitoxin", "addiction module",
                                   "killer protein", "antitoxin",
                                   "type ii toxin", "plasmid stabilization"))


def _is_virulence_gene(gene_name: str, product: str) -> bool:
    """Check if a gene is a virulence factor."""
    vir_prefixes = ("vir", "iro", "iuc", "iut", "ybt", "fyua", "sit", "kfu",
                    "clb", "cnf", "hly", "tsh", "iss", "cia", "cva", "cvab",
                    "etsa", "etsb", "pic", "sat", "pet", "esp", "eae",
                    "afa", "pap", "fim", "sfa", "iha", "hra")
    g = gene_name.lower()
    if any(g.startswith(p) for p in vir_prefixes):
        return True
    p = product.lower()
    return any(kw in p for kw in ("virulence", "hemolysin", "siderophore",
                                   "aerobactin", "yersiniabactin",
                                   "colicin", "invasin", "adhesin"))


def _is_qac_gene(gene_name: str, product: str) -> bool:
    """Check if a gene confers quaternary ammonium compound resistance."""
    qac_names = {"qace", "qaca", "qacb", "qacc", "qacd", "qacf", "qacg",
                 "qach", "qacj", "qacr", "smr", "emre", "sugE"}
    g = gene_name.lower()
    if g in qac_names or g.startswith("qac"):
        return True
    p = product.lower()
    return any(kw in p for kw in ("quaternary ammonium", "qac efflux",
                                   "small multidrug resistance"))


def _is_metal_gene(gene_name: str, product: str) -> bool:
    """Check if a gene confers heavy metal resistance."""
    metal_prefixes = ("mer", "ars", "pco", "sil", "ter", "cop", "cus", "czc",
                      "cad", "pbr", "chr", "nik")
    g = gene_name.lower()
    if any(g.startswith(p) for p in metal_prefixes):
        return True
    p = product.lower()
    return any(kw in p for kw in ("mercury", "arsenic", "copper", "silver",
                                   "tellurium", "lead", "cadmium", "chromate",
                                   "heavy metal"))


def classify_gene(gene_name: str, product: str) -> str:
    """Classify a gene into a functional category. Returns the category name."""
    if _is_amr_gene(gene_name, product):
        return "AMR"
    if _is_mob_gene(gene_name, product):
        return "MOB"
    if _is_ta_gene(gene_name, product):
        return "TA"
    if _is_virulence_gene(gene_name, product):
        return "VIR"
    if _is_qac_gene(gene_name, product):
        return "QAC"
    if _is_metal_gene(gene_name, product):
        return "Metal"
    return "CDS"


def _nice_tick_interval(length):
    """Calculate a nice tick interval for the backbone labels."""
    if length <= 5000:
        return 1000
    elif length <= 20000:
        return 5000
    elif length <= 100000:
        return 10000
    elif length <= 500000:
        return 50000
    else:
        return 100000


def _format_size(bp):
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f} Mb"
    elif bp >= 1000:
        return f"{bp / 1000:.1f} kb"
    return f"{bp} bp"


# ---------------------------------------------------------------------------
# Inc Groups & Mobility Charts
# ---------------------------------------------------------------------------

def make_inc_group_chart():
    """Bar chart of incompatibility group (replicon type) distribution."""
    rep = mob_typing.get("rep_types", {})
    if not rep:
        return go.Figure()
    df = pd.DataFrame(
        sorted(rep.items(), key=lambda x: x[1], reverse=True)[:20],
        columns=["Inc Group / Replicon Type", "Count"],
    )
    fig = px.bar(
        df, x="Count", y="Inc Group / Replicon Type", orientation="h",
        color="Count", color_continuous_scale="Tealgrn",
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=500,
        yaxis=dict(autorange="reversed"), coloraxis_showscale=False,
        xaxis_title="Number of Plasmids", yaxis_title="",
    )
    return fig


def make_mobility_chart():
    """Donut chart of predicted mobility distribution."""
    mob = mob_typing.get("mobility", {})
    if not mob:
        return go.Figure()
    color_map = {
        "conjugative": COLORS["accent3"],
        "mobilizable": COLORS["accent"],
        "non-mobilizable": COLORS["accent4"],
    }
    labels = list(mob.keys())
    values = list(mob.values())
    colors = [color_map.get(l, "#64748b") for l in labels]

    fig = go.Figure(go.Pie(
        labels=[l.capitalize() for l in labels], values=values,
        hole=0.55, marker_colors=colors,
        textinfo="label+percent", textfont_size=13,
        hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>Percent: %{percent}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        showlegend=False,
    )
    return fig


def make_relaxase_chart():
    """Bar chart of relaxase type distribution."""
    relax = mob_typing.get("relaxase_types", {})
    if not relax:
        return go.Figure()
    df = pd.DataFrame(
        sorted(relax.items(), key=lambda x: x[1], reverse=True),
        columns=["Relaxase Type", "Count"],
    )
    fig = px.bar(
        df, x="Relaxase Type", y="Count",
        color="Relaxase Type",
        color_discrete_sequence=px.colors.qualitative.Vivid,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis_title="", yaxis_title="Number of Plasmids",
        showlegend=False,
    )
    return fig


def make_mpf_chart():
    """Bar chart of MPF (Mating Pair Formation) type distribution."""
    mpf = mob_typing.get("mpf_types", {})
    if not mpf:
        return go.Figure()
    df = pd.DataFrame(
        sorted(mpf.items(), key=lambda x: x[1], reverse=True),
        columns=["MPF Type", "Count"],
    )
    fig = px.bar(
        df, x="MPF Type", y="Count",
        color="MPF Type",
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis_title="", yaxis_title="Number of Plasmids",
        showlegend=False,
    )
    return fig


# ---------------------------------------------------------------------------
# Correlation Charts
# ---------------------------------------------------------------------------

def make_amr_inc_heatmap():
    """Heatmap of AMR drug class vs Inc group co-occurrence."""
    amr_inc = correlations.get("amr_vs_inc", {})
    if not amr_inc:
        return go.Figure()

    all_inc = {}
    for dc_data in amr_inc.values():
        for inc, cnt in dc_data.items():
            all_inc[inc] = all_inc.get(inc, 0) + cnt
    top_inc = sorted(all_inc, key=all_inc.get, reverse=True)[:12]
    drug_classes = sorted(amr_inc.keys(),
                          key=lambda dc: sum(amr_inc[dc].values()), reverse=True)[:12]

    z = [[amr_inc[dc].get(inc, 0) for inc in top_inc] for dc in drug_classes]

    fig = go.Figure(go.Heatmap(
        z=z, x=top_inc, y=drug_classes,
        colorscale="YlOrRd", texttemplate="%{z}",
        textfont={"size": 10},
        hovertemplate="Drug: <b>%{y}</b><br>Inc Group: <b>%{x}</b><br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=480,
        xaxis=dict(tickangle=-45, side="bottom"),
        yaxis=dict(autorange="reversed"),
    )
    return fig


def make_amr_mobility_chart():
    """Grouped bar chart of AMR drug class vs predicted mobility."""
    amr_mob = correlations.get("amr_vs_mobility", {})
    if not amr_mob:
        return go.Figure()

    drug_classes = sorted(amr_mob.keys(),
                          key=lambda dc: sum(amr_mob[dc].values()), reverse=True)[:12]
    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]

    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        values = [amr_mob[dc].get(mob, 0) for dc in drug_classes]
        fig.add_trace(go.Bar(
            x=drug_classes, y=values, name=mob.capitalize(),
            marker_color=color,
            hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>",
        ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=420,
        barmode="group", xaxis_tickangle=-45,
        xaxis_title="", yaxis_title="Number of Plasmids",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_heavy_metal_genes_chart():
    """Bar chart of heavy metal resistance gene frequencies."""
    hm = correlations.get("heavy_metal_genes", {})
    if not hm:
        return go.Figure()

    metal_map = {}
    for gene in sorted(hm.keys()):
        g = gene.lower()
        if g.startswith("mer"): metal_map[gene] = "Mercury"
        elif g.startswith("ars"): metal_map[gene] = "Arsenic"
        elif g.startswith("pco"): metal_map[gene] = "Copper"
        elif g.startswith("sil"): metal_map[gene] = "Silver"
        elif g.startswith("ter"): metal_map[gene] = "Tellurium"
        else: metal_map[gene] = "Other"

    df = pd.DataFrame([
        {"Gene": g, "Count": c, "Metal": metal_map.get(g, "Other")}
        for g, c in sorted(hm.items(), key=lambda x: x[1], reverse=True)
    ])

    metal_colors = {
        "Mercury": "#ef4444", "Arsenic": "#f59e0b", "Copper": "#f97316",
        "Silver": "#94a3b8", "Tellurium": "#8b5cf6", "Other": "#64748b",
    }
    fig = px.bar(df, x="Gene", y="Count", color="Metal",
                 color_discrete_map=metal_colors)
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis_tickangle=-45, xaxis_title="", yaxis_title="Number of Plasmids",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_heavy_metal_inc_heatmap():
    """Heatmap of heavy metal resistance vs Inc group."""
    hm_inc = correlations.get("heavy_metal_vs_inc", {})
    if not hm_inc:
        return go.Figure()

    metals = sorted(hm_inc.keys())
    all_inc = {}
    for m_data in hm_inc.values():
        for inc, cnt in m_data.items():
            all_inc[inc] = all_inc.get(inc, 0) + cnt
    top_inc = sorted(all_inc, key=all_inc.get, reverse=True)[:10]
    z = [[hm_inc[m].get(inc, 0) for inc in top_inc] for m in metals]

    fig = go.Figure(go.Heatmap(
        z=z, x=top_inc, y=metals,
        colorscale="Purples", texttemplate="%{z}", textfont={"size": 11},
        hovertemplate="Metal: <b>%{y}</b><br>Inc: <b>%{x}</b><br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


def make_heavy_metal_mobility_chart():
    """Stacked bar: heavy metal resistance by mobility type."""
    hm_mob = correlations.get("heavy_metal_vs_mobility", {})
    if not hm_mob:
        return go.Figure()

    metals = sorted(hm_mob.keys(), key=lambda m: sum(hm_mob[m].values()), reverse=True)
    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]

    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        vals = [hm_mob[m].get(mob, 0) for m in metals]
        fig.add_trace(go.Bar(x=metals, y=vals, name=mob.capitalize(), marker_color=color))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        barmode="stack", xaxis_title="", yaxis_title="Plasmid Count",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_phage_distribution_chart():
    """Side-by-side bars for phage elements and transposase counts."""
    phage = correlations.get("phage_distribution", {})
    trans = correlations.get("transposase_distribution", {})
    if not phage and not trans:
        return go.Figure()

    fig = go.Figure()
    if phage:
        fig.add_trace(go.Bar(
            x=list(phage.keys()), y=list(phage.values()),
            name="Phage / prophage elements", marker_color=COLORS["accent5"],
        ))
    if trans:
        fig.add_trace(go.Bar(
            x=list(trans.keys()), y=list(trans.values()),
            name="Transposases", marker_color=COLORS["accent2"],
        ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        barmode="group", xaxis_title="Elements per plasmid",
        yaxis_title="Number of Plasmids",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_phage_mobility_chart():
    """Bar chart: mean phage/integrase elements by mobility class."""
    phage_mob = correlations.get("phage_vs_mobility", {})
    if not phage_mob:
        return go.Figure()

    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]
    means, counts, labels, colors = [], [], [], []
    for mob, color in zip(mob_types, mob_colors):
        if mob in phage_mob:
            labels.append(mob.capitalize())
            means.append(phage_mob[mob]["mean_phage_elements"])
            counts.append(phage_mob[mob]["count"])
            colors.append(color)

    fig = go.Figure(go.Bar(
        x=labels, y=means, marker_color=colors,
        text=[f"{m:.1f}" for m in means], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Mean elements: %{y:.1f}<br>n = %{customdata}<extra></extra>",
        customdata=counts,
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=30, b=10, l=10, r=10), height=350,
        xaxis_title="", yaxis_title="Mean phage + integrase elements",
    )
    return fig


def make_pmlst_scheme_chart():
    """Bar chart of pMLST scheme distribution."""
    schemes = correlations.get("pmlst_schemes", {})
    if not schemes:
        return go.Figure()

    df = pd.DataFrame(
        sorted(schemes.items(), key=lambda x: x[1], reverse=True),
        columns=["pMLST Scheme", "Count"],
    )
    fig = px.bar(df, x="pMLST Scheme", y="Count",
                 color="Count", color_continuous_scale="Tealgrn")
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis_title="", yaxis_title="Number of Plasmids",
        coloraxis_showscale=False,
    )
    return fig


def make_pmlst_alleles_chart():
    """Bar chart of pMLST allele frequencies."""
    alleles = correlations.get("pmlst_alleles", {})
    if not alleles:
        return go.Figure()

    df = pd.DataFrame(
        sorted(alleles.items(), key=lambda x: x[1], reverse=True)[:20],
        columns=["Allele", "Count"],
    )
    fig = px.bar(df, x="Allele", y="Count",
                 color_discrete_sequence=[COLORS["accent"]])
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis_title="", yaxis_title="Number of Plasmids",
        xaxis_tickangle=-45,
    )
    return fig


def make_pmlst_mobility_chart():
    """Stacked bar: pMLST scheme vs mobility."""
    pmlst_mob = correlations.get("pmlst_vs_mobility", {})
    if not pmlst_mob:
        return go.Figure()

    schemes = sorted(pmlst_mob.keys(),
                     key=lambda s: sum(pmlst_mob[s].values()), reverse=True)
    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]

    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        vals = [pmlst_mob[s].get(mob, 0) for s in schemes]
        fig.add_trace(go.Bar(x=schemes, y=vals, name=mob.capitalize(), marker_color=color))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        barmode="stack", xaxis_title="", yaxis_title="Plasmid Count",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


# ---------------------------------------------------------------------------
# Virulence, TA, QAC correlation charts
# ---------------------------------------------------------------------------

def _make_feature_inc_chart(data_key, title_label, bar_color):
    """Bar chart: Inc groups carrying a specific feature category."""
    inc_data = correlations.get(data_key, {})
    if not inc_data:
        return go.Figure()
    df = pd.DataFrame(
        sorted(inc_data.items(), key=lambda x: x[1], reverse=True)[:15],
        columns=["Inc Group", "Plasmids"],
    )
    fig = px.bar(df, x="Plasmids", y="Inc Group", orientation="h",
                 color_discrete_sequence=[bar_color])
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=420,
        yaxis=dict(autorange="reversed"),
        xaxis_title=f"Plasmids with {title_label}", yaxis_title="",
    )
    return fig


def _make_feature_mobility_chart(data_key, title_label):
    """Donut: mobility distribution of plasmids carrying a feature."""
    mob_data = correlations.get(data_key, {})
    if not mob_data or sum(mob_data.values()) == 0:
        return go.Figure()
    labels = [k.capitalize() for k in mob_data.keys()]
    values = list(mob_data.values())
    colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]
    fig = go.Figure(go.Pie(
        labels=labels, values=values, hole=0.55,
        marker_colors=colors[:len(labels)],
        textinfo="label+percent", textfont_size=12,
        hovertemplate="<b>%{label}</b><br>Count: %{value:,}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350, showlegend=False,
    )
    return fig


def _make_gene_frequency_chart(data_key, title_label, bar_color):
    """Bar chart: top gene names for a feature category."""
    gene_data = correlations.get(data_key, {})
    if not gene_data:
        return go.Figure()
    df = pd.DataFrame(
        sorted(gene_data.items(), key=lambda x: x[1], reverse=True)[:20],
        columns=["Gene", "Count"],
    )
    fig = px.bar(df, x="Gene", y="Count", color_discrete_sequence=[bar_color])
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350,
        xaxis_tickangle=-45, xaxis_title="", yaxis_title="Plasmid Count",
    )
    return fig


# ---------------------------------------------------------------------------
# Geography charts
# ---------------------------------------------------------------------------

def make_global_map():
    """Scatter map of plasmid reports coloured by predicted mobility."""
    geo = db.geo_plasmid_data()
    if not geo:
        return go.Figure()

    df = pd.DataFrame(geo)
    mob_colors = {
        "conjugative": COLORS["accent3"],
        "mobilizable": COLORS["accent"],
        "non-mobilizable": COLORS["accent4"],
    }
    fig = px.scatter_geo(
        df, lat="lat", lon="lng", color="predicted_mobility",
        hover_name="NUCCORE_ACC",
        hover_data={"country": True, "rep_type": True, "predicted_mobility": True,
                    "lat": False, "lng": False},
        color_discrete_map=mob_colors,
        opacity=0.5, size_max=8,
        projection="natural earth",
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=550,
        geo=dict(
            showland=True, landcolor="#f0ebe3",
            showocean=True, oceancolor="#e8f4f8",
            showcountries=True, countrycolor="#cbd5e1",
            showcoastlines=True, coastlinecolor="#94a3b8",
        ),
        legend=dict(orientation="h", yanchor="top", y=-0.02,
                    xanchor="center", x=0.5, font_size=11),
    )
    return fig


def make_country_mobility_chart():
    """Stacked bar: mobility distribution by top 20 countries."""
    geo_summary = db.geo_country_summary()
    if not geo_summary:
        return go.Figure()

    # Aggregate by country
    country_mob = {}
    for r in geo_summary:
        c = r["country"]
        m = r["predicted_mobility"] or "unknown"
        country_mob.setdefault(c, {})
        country_mob[c][m] = country_mob[c].get(m, 0) + r["cnt"]

    # Top 20 by total
    top = sorted(country_mob.items(),
                 key=lambda x: sum(x[1].values()), reverse=True)[:20]
    countries = [c for c, _ in top]

    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]

    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        vals = [country_mob[c].get(mob, 0) for c in countries]
        fig.add_trace(go.Bar(x=countries, y=vals, name=mob.capitalize(),
                             marker_color=color))
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        barmode="stack", xaxis_title="", yaxis_title="Plasmid Count",
        xaxis_tickangle=-45,
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_country_inc_chart():
    """Heatmap: top Inc groups vs top countries."""
    country_inc = db.geo_country_inc_top(15)
    if not country_inc:
        return go.Figure()

    # Get top Inc groups across all countries
    all_inc = {}
    for country, incs in country_inc.items():
        for inc, cnt in incs.items():
            all_inc[inc] = all_inc.get(inc, 0) + cnt
    top_inc = sorted(all_inc, key=all_inc.get, reverse=True)[:12]
    countries = sorted(country_inc.keys(),
                       key=lambda c: sum(country_inc[c].values()), reverse=True)

    z = [[country_inc[c].get(inc, 0) for inc in top_inc] for c in countries]

    fig = go.Figure(go.Heatmap(
        z=z, x=top_inc, y=countries,
        colorscale="YlGnBu", texttemplate="%{z}",
        textfont={"size": 9},
        hovertemplate="Country: <b>%{y}</b><br>Inc: <b>%{x}</b><br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=500,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


def make_temporal_animation():
    """Animated scatter map: plasmid reports by year."""
    temporal = db.geo_temporal_data()
    if not temporal:
        return go.Figure()

    df = pd.DataFrame(temporal)
    df["year"] = df["year"].astype(str)

    # Cumulative counts per country per year
    all_years = sorted(df["year"].unique())
    cumulative_rows = []
    country_totals = {}
    for year in all_years:
        year_data = df[df["year"] == year]
        for _, row in year_data.iterrows():
            c = row["country"]
            country_totals[c] = country_totals.get(c, 0) + row["cnt"]
        # Add all countries seen so far for this frame
        for c, total in country_totals.items():
            lat_row = df[df["country"] == c].iloc[0]
            cumulative_rows.append({
                "year": year, "country": c,
                "lat": lat_row["lat"], "lng": lat_row["lng"],
                "cumulative": total,
            })

    cum_df = pd.DataFrame(cumulative_rows)

    fig = px.scatter_geo(
        cum_df, lat="lat", lon="lng",
        size="cumulative", color="cumulative",
        hover_name="country",
        hover_data={"cumulative": True, "year": True, "lat": False, "lng": False},
        animation_frame="year",
        color_continuous_scale="YlOrRd",
        size_max=40,
        projection="natural earth",
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=30, b=10, l=10, r=10), height=600,
        geo=dict(
            showland=True, landcolor="#f0ebe3",
            showocean=True, oceancolor="#e8f4f8",
            showcountries=True, countrycolor="#cbd5e1",
            showcoastlines=True, coastlinecolor="#94a3b8",
        ),
        coloraxis_colorbar=dict(title="Plasmids"),
    )
    fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 800
    fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 400
    return fig


def _build_animated_map(raw_data, title_suffix="", color_scale="YlOrRd"):
    """Generic cumulative animated scatter map builder."""
    if not raw_data:
        return go.Figure()
    df = pd.DataFrame(raw_data)
    df["year"] = df["year"].astype(str)

    all_years = sorted(df["year"].unique())
    cumulative_rows = []
    country_totals = {}
    for year in all_years:
        year_data = df[df["year"] == year]
        for _, row in year_data.iterrows():
            c = row["country"]
            country_totals[c] = country_totals.get(c, 0) + row["cnt"]
        for c, total in country_totals.items():
            lat_row = df[df["country"] == c].iloc[0]
            cumulative_rows.append({
                "year": year, "country": c,
                "lat": lat_row["lat"], "lng": lat_row["lng"],
                "cumulative": total,
            })

    cum_df = pd.DataFrame(cumulative_rows)
    fig = px.scatter_geo(
        cum_df, lat="lat", lon="lng",
        size="cumulative", color="cumulative",
        hover_name="country",
        hover_data={"cumulative": True, "year": True, "lat": False, "lng": False},
        animation_frame="year",
        color_continuous_scale=color_scale,
        size_max=40, projection="natural earth",
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color=COLORS["text"],
        margin=dict(t=30, b=10, l=10, r=10), height=550,
        geo=dict(showland=True, landcolor="#f0ebe3",
                 showocean=True, oceancolor="#e8f4f8",
                 showcountries=True, countrycolor="#cbd5e1",
                 showcoastlines=True, coastlinecolor="#94a3b8"),
        coloraxis_colorbar=dict(title="Plasmids"),
    )
    if fig.layout.updatemenus:
        fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 800
        fig.layout.updatemenus[0].buttons[0].args[1]["transition"]["duration"] = 400
    return fig


def make_host_distribution_chart():
    """Donut chart: plasmid host/source categories."""
    host_counts = db.host_category_counts()
    if not host_counts:
        return go.Figure()
    host_colors = {
        "Human": "#e74c3c", "Animal": "#f39c12", "Soil": "#8b4513",
        "Water": "#3498db", "Food": "#2ecc71", "Environment": "#9b59b6",
    }
    labels = list(host_counts.keys())
    values = list(host_counts.values())
    colors = [host_colors.get(l, "#95a5a6") for l in labels]
    fig = go.Figure(go.Pie(
        labels=labels, values=values, hole=0.5,
        marker_colors=colors, textinfo="label+percent", textfont_size=12,
        hovertemplate="<b>%{label}</b><br>Plasmids: %{value:,}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400, showlegend=False,
    )
    return fig


def make_host_mobility_chart():
    """Stacked bar: mobility by host category."""
    host_mob = db.host_vs_mobility()
    if not host_mob:
        return go.Figure()
    hosts = sorted(host_mob.keys(), key=lambda h: sum(host_mob[h].values()), reverse=True)
    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]
    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        fig.add_trace(go.Bar(
            x=hosts, y=[host_mob[h].get(mob, 0) for h in hosts],
            name=mob.capitalize(), marker_color=color,
        ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        barmode="stack", xaxis_title="", yaxis_title="Plasmid Count",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_host_inc_heatmap():
    """Heatmap: top Inc groups vs host categories."""
    host_inc = db.host_vs_inc()
    if not host_inc:
        return go.Figure()
    all_inc = {}
    for cat_data in host_inc.values():
        for inc, cnt in cat_data.items():
            all_inc[inc] = all_inc.get(inc, 0) + cnt
    top_inc = sorted(all_inc, key=all_inc.get, reverse=True)[:12]
    hosts = sorted(host_inc.keys(), key=lambda h: sum(host_inc[h].values()), reverse=True)
    z = [[host_inc[h].get(inc, 0) for inc in top_inc] for h in hosts]
    fig = go.Figure(go.Heatmap(
        z=z, x=top_inc, y=hosts,
        colorscale="YlGnBu", texttemplate="%{z}", textfont={"size": 9},
        hovertemplate="Host: <b>%{y}</b><br>Inc: <b>%{x}</b><br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


def make_host_amr_heatmap():
    """Heatmap: top AMR drug classes vs host categories."""
    host_amr = db.host_vs_amr_class()
    if not host_amr:
        return go.Figure()
    # Get top drug classes
    all_dc = {}
    skip = {"", "COPPER", "COPPER/SILVER", "MERCURY", "SILVER",
            "TELLURIUM", "ARSENIC", "ZINC", "LEAD", "CADMIUM", "CHROMIUM"}
    for cat_data in host_amr.values():
        for dc, cnt in cat_data.items():
            if dc.upper() not in skip:
                all_dc[dc] = all_dc.get(dc, 0) + cnt
    top_dc = sorted(all_dc, key=all_dc.get, reverse=True)[:12]
    hosts = sorted(host_amr.keys(), key=lambda h: sum(host_amr[h].values()), reverse=True)
    z = [[host_amr[h].get(dc, 0) for dc in top_dc] for h in hosts]
    fig = go.Figure(go.Heatmap(
        z=z, x=top_dc, y=hosts,
        colorscale="YlOrRd", texttemplate="%{z}", textfont={"size": 9},
        hovertemplate="Host: <b>%{y}</b><br>Drug: <b>%{x}</b><br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


# ---------------------------------------------------------------------------
# Analytics charts
# ---------------------------------------------------------------------------

def make_matched_comparison():
    """Heatmap: conjugative % per species x country (matched comparison)."""
    rows = db.analytics_matched_comparison()
    if not rows:
        return go.Figure()

    from collections import defaultdict
    data = defaultdict(lambda: defaultdict(lambda: {"conj": 0, "total": 0}))
    for r in rows:
        d = data[r["genus"]][r["country"]]
        d["total"] += r["cnt"]
        if r["mobility"] == "conjugative":
            d["conj"] += r["cnt"]

    genera = sorted(data.keys())
    countries = sorted({c for g in data.values() for c in g.keys()},
                       key=lambda c: -sum(data[g][c]["total"] for g in genera))

    z = []
    text = []
    for genus in genera:
        row_z = []
        row_t = []
        for country in countries:
            d = data[genus][country]
            if d["total"] >= 10:
                pct = d["conj"] / d["total"] * 100
                row_z.append(round(pct, 1))
                row_t.append(f"{pct:.0f}%<br>n={d['total']}")
            else:
                row_z.append(None)
                row_t.append("")
        z.append(row_z)
        text.append(row_t)

    fig = go.Figure(go.Heatmap(
        z=z, x=countries, y=genera,
        colorscale="RdYlGn", zmid=50,
        text=text, texttemplate="%{text}", textfont={"size": 9},
        hovertemplate="<b>%{y}</b> in %{x}<br>Conjugative: %{z:.1f}%<extra></extra>",
        colorbar=dict(title="% Conjugative"),
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


def make_rarefaction_chart():
    """Rarefaction curves: conjugative fraction vs sample size per country."""
    results = db.analytics_rarefaction(n_iter=50)
    if not results:
        return go.Figure()

    df = pd.DataFrame(results)
    fig = go.Figure()
    colors = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel
    for i, country in enumerate(df["country"].unique()):
        cdf = df[df["country"] == country].sort_values("n")
        fig.add_trace(go.Scatter(
            x=cdf["n"], y=cdf["conj_mean"],
            mode="lines+markers", name=country,
            line=dict(color=colors[i % len(colors)]),
            error_y=dict(type="data", array=cdf["conj_std"].tolist(),
                         visible=True, thickness=1),
            hovertemplate=(f"<b>{country}</b><br>"
                           "n=%{x}<br>Conjugative: %{y:.1%}<extra></extra>"),
        ))

    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        xaxis_title="Sample size (n)", yaxis_title="Conjugative fraction",
        yaxis=dict(tickformat=".0%"),
        legend=dict(font_size=10),
    )
    return fig


def make_feature_importance_chart():
    """Random Forest feature importance for predicting mobility."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import LabelEncoder
    import numpy as np

    rows = db.analytics_feature_matrix()
    if not rows or len(rows) < 100:
        return go.Figure()

    df = pd.DataFrame(rows)

    # Encode categoricals
    features = ["country", "genus", "host_category", "year", "length", "gc"]
    encoders = {}
    X_data = {}
    for col in ["country", "genus", "host_category", "year"]:
        le = LabelEncoder()
        X_data[col] = le.fit_transform(df[col].fillna("Unknown"))
        encoders[col] = le
    X_data["length"] = df["length"].fillna(0).values.astype(float)
    X_data["gc"] = df["gc"].fillna(0).values.astype(float)

    # Add Inc group features (top 10)
    top_inc = ["IncFIB", "IncFII", "IncFIA", "IncI1", "IncN",
               "IncHI2", "IncC", "IncX4", "ColRNAI", "IncR"]
    for inc in top_inc:
        X_data[f"has_{inc}"] = df["rep_type"].fillna("").str.contains(inc).astype(int).values
        features.append(f"has_{inc}")

    X = pd.DataFrame(X_data)[features]
    y = LabelEncoder().fit_transform(df["mobility"])

    # Train Random Forest
    rf = RandomForestClassifier(n_estimators=100, max_depth=10,
                                random_state=42, n_jobs=-1)
    rf.fit(X, y)
    importances = rf.feature_importances_

    # Clean feature names for display
    display_names = []
    for f in features:
        if f.startswith("has_"):
            display_names.append(f[4:])
        else:
            display_names.append(f.capitalize())

    imp_df = pd.DataFrame({
        "Feature": display_names,
        "Importance": importances,
    }).sort_values("Importance", ascending=True)

    fig = px.bar(imp_df, x="Importance", y="Feature", orientation="h",
                 color="Importance", color_continuous_scale="Viridis")
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=500,
        coloraxis_showscale=False,
        xaxis_title="Feature Importance (Gini)", yaxis_title="",
    )
    return fig


def make_simpson_paradox_table():
    """Table showing Simpson's paradox and high species-divergence cases."""
    paradoxes = db.analytics_simpson_paradox()
    if not paradoxes:
        return html.P("No significant species-level divergence detected. "
                      "Regional mobility differences are consistent across species.",
                      className="text-muted")

    rows = []
    for p in paradoxes:
        species_detail = ", ".join(
            f"{sp}: {pct}%" for sp, pct in list(p["species"].items())[:4]
        )
        flag = "PARADOX" if p.get("paradox") else f"{p.get('spread', 0)}pp spread"
        rows.append(html.Tr([
            html.Td(p["inc_group"]),
            html.Td(f"{p['overall_conj_pct']}%"),
            html.Td(f"{p['total']:,}"),
            html.Td(flag, style={"color": "#e74c3c" if p.get("paradox") else "#f39c12",
                                 "fontWeight": "600"}),
            html.Td(species_detail, style={"fontSize": "0.85rem"}),
        ]))

    return html.Table(className="analytics-table", children=[
        html.Thead(html.Tr([
            html.Th("Inc Group"),
            html.Th("Overall % Conj."),
            html.Th("n"),
            html.Th("Status"),
            html.Th("Per-species % Conjugative"),
        ])),
        html.Tbody(rows),
    ])


def make_umap_chart(color_by="mobility"):
    """UMAP scatter plot of plasmid feature space."""
    umap_path = os.path.join(os.path.dirname(__file__), "data", "umap_embedding.npz")
    if not os.path.exists(umap_path):
        return go.Figure()

    import numpy as np
    data = np.load(umap_path, allow_pickle=True)
    df = pd.DataFrame({
        "UMAP1": data["x"], "UMAP2": data["y"],
        "Mobility": data["mobility"], "Genus": data["genus"],
        "Country": data["country"], "Host": data["host"],
        "Accession": data["acc"],
    })

    color_map = {
        "mobility": {"col": "Mobility", "palette": {
            "conjugative": COLORS["accent3"],
            "mobilizable": COLORS["accent"],
            "non-mobilizable": COLORS["accent4"],
        }},
        "genus": {"col": "Genus", "palette": None},
        "host": {"col": "Host", "palette": None},
    }

    cfg = color_map.get(color_by, color_map["mobility"])
    col = cfg["col"]

    if cfg["palette"]:
        fig = px.scatter(df, x="UMAP1", y="UMAP2", color=col,
                         color_discrete_map=cfg["palette"],
                         hover_data=["Accession", "Genus", "Mobility"],
                         opacity=0.5)
    else:
        # Limit to top 8 categories
        top = df[col].value_counts().head(8).index.tolist()
        df_plot = df[df[col].isin(top)]
        fig = px.scatter(df_plot, x="UMAP1", y="UMAP2", color=col,
                         color_discrete_sequence=px.colors.qualitative.Set2,
                         hover_data=["Accession", "Mobility"],
                         opacity=0.5)

    fig.update_traces(marker_size=3)
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=550,
        xaxis_title="UMAP 1", yaxis_title="UMAP 2",
        legend=dict(font_size=10),
    )
    return fig


def make_xgboost_shap_chart():
    """XGBoost + SHAP: feature importance with interaction effects."""
    import xgboost as xgb
    import shap
    import numpy as np
    from sklearn.preprocessing import LabelEncoder

    rows = db.analytics_feature_matrix()
    if not rows or len(rows) < 500:
        return go.Figure(), None

    df = pd.DataFrame(rows)

    # Encode features
    features = []
    X_cols = {}
    for col in ["country", "genus", "host_category", "year"]:
        le = LabelEncoder()
        X_cols[col] = le.fit_transform(df[col].fillna("Unknown"))
        features.append(col)

    X_cols["length"] = np.log1p(df["length"].fillna(0).values.astype(float))
    features.append("length")
    X_cols["gc"] = df["gc"].fillna(0).values.astype(float)
    features.append("gc")

    # Top Inc groups as binary features
    top_inc = ["IncFIB", "IncFII", "IncFIA", "IncI1", "IncN",
               "IncHI2", "IncC", "IncX4", "ColRNAI", "IncR"]
    for inc in top_inc:
        X_cols[inc] = df["rep_type"].fillna("").str.contains(inc).astype(int).values
        features.append(inc)

    X = pd.DataFrame(X_cols)[features]
    y_le = LabelEncoder()
    y = y_le.fit_transform(df["mobility"])

    # Train XGBoost
    model = xgb.XGBClassifier(
        n_estimators=200, max_depth=6, learning_rate=0.1,
        subsample=0.8, colsample_bytree=0.8,
        random_state=42, n_jobs=-1, eval_metric="mlogloss",
    )
    model.fit(X, y)

    # SHAP values (use TreeExplainer for speed)
    explainer = shap.TreeExplainer(model)
    # Sample for speed
    X_sample = X.sample(min(5000, len(X)), random_state=42)
    shap_values = explainer.shap_values(X_sample)

    # Mean absolute SHAP per feature (across all classes)
    sv = np.array(shap_values)
    if sv.ndim == 3:
        # Shape (n_samples, n_features, n_classes) — average over samples and classes
        mean_shap = np.mean(np.abs(sv), axis=(0, 2))
    elif sv.ndim == 2:
        mean_shap = np.abs(sv).mean(axis=0)
    else:
        mean_shap = np.abs(sv).flatten()
    mean_shap = mean_shap[:len(features)]  # safety trim

    # Display names
    display = [f.replace("_", " ").capitalize() for f in features]

    imp_df = pd.DataFrame({
        "Feature": display, "SHAP Importance": mean_shap,
    }).sort_values("SHAP Importance", ascending=True)

    fig = px.bar(imp_df, x="SHAP Importance", y="Feature", orientation="h",
                 color="SHAP Importance", color_continuous_scale="Viridis")
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=500,
        coloraxis_showscale=False,
        xaxis_title="Mean |SHAP value|", yaxis_title="",
    )

    # Accuracy info
    from sklearn.metrics import accuracy_score
    acc = accuracy_score(y, model.predict(X))

    return fig, round(acc * 100, 1)


def make_temporal_trends_chart():
    """Line chart: AMR drug class prevalence trends by year."""
    rows = db.q("""
        SELECT SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
               a.drug_class, COUNT(DISTINCT n.NUCCORE_ACC) AS cnt
        FROM nuccore n
        JOIN amr a ON n.NUCCORE_ACC = a.NUCCORE_ACC
        WHERE a.drug_class IN ('BETA-LACTAM', 'AMINOGLYCOSIDE', 'TETRACYCLINE',
                               'SULFONAMIDE', 'QUINOLONE', 'TRIMETHOPRIM',
                               'PHENICOL', 'MERCURY', 'QUATERNARY AMMONIUM')
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2010
        GROUP BY year, a.drug_class
        ORDER BY year
    """)
    if not rows:
        return go.Figure()

    # Normalize by total plasmids per year
    year_totals = db.q("""
        SELECT SUBSTR(NUCCORE_CreateDate, 1, 4) AS year, COUNT(*) AS total
        FROM nuccore
        WHERE LENGTH(NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2010
        GROUP BY year
    """)
    totals = {r["year"]: r["total"] for r in year_totals}

    df = pd.DataFrame(rows)
    df["prevalence"] = df.apply(
        lambda r: r["cnt"] / totals.get(r["year"], 1) * 100, axis=1)

    fig = px.line(
        df, x="year", y="prevalence", color="drug_class",
        markers=True,
        color_discrete_sequence=px.colors.qualitative.Set2 + px.colors.qualitative.Pastel,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        xaxis_title="Year", yaxis_title="% of plasmids carrying resistance",
        legend=dict(title="Drug Class", font_size=9),
    )
    return fig


def make_inc_trends_chart():
    """Line chart: Inc group prevalence trends by year."""
    rows = db.q("""
        SELECT SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year, t.rep_type
        FROM nuccore n
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE t.rep_type != ''
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2010
    """)
    if not rows:
        return go.Figure()

    year_totals = db.q("""
        SELECT SUBSTR(NUCCORE_CreateDate, 1, 4) AS year, COUNT(*) AS total
        FROM nuccore WHERE LENGTH(NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2010
        GROUP BY year
    """)
    totals = {r["year"]: r["total"] for r in year_totals}

    top_inc = ["IncFIB", "IncFII", "IncFIA", "IncI1", "IncN", "IncHI2", "ColRNAI", "IncC"]
    year_inc = {}
    for r in rows:
        year = r["year"]
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt in top_inc:
                year_inc.setdefault(year, {}).setdefault(rt, 0)
                year_inc[year][rt] += 1

    plot_rows = []
    for year in sorted(year_inc.keys()):
        total = totals.get(year, 1)
        for inc in top_inc:
            cnt = year_inc[year].get(inc, 0)
            plot_rows.append({"year": year, "Inc Group": inc,
                              "prevalence": cnt / total * 100})

    df = pd.DataFrame(plot_rows)
    fig = px.line(df, x="year", y="prevalence", color="Inc Group", markers=True,
                  color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        xaxis_title="Year", yaxis_title="% of plasmids",
        legend=dict(title="Inc Group", font_size=9),
    )
    return fig


def make_amr_cooccurrence_chart():
    """Heatmap of AMR gene co-occurrence on the same plasmid."""
    rows = db.amr_cooccurrence(min_count=100)
    if not rows:
        return go.Figure()

    # Build matrix from top gene pairs
    genes = set()
    for r in rows[:50]:
        genes.add(r["gene1"])
        genes.add(r["gene2"])
    genes = sorted(genes)

    matrix = {g: {g2: 0 for g2 in genes} for g in genes}
    for r in rows:
        if r["gene1"] in matrix and r["gene2"] in matrix:
            matrix[r["gene1"]][r["gene2"]] = r["cnt"]
            matrix[r["gene2"]][r["gene1"]] = r["cnt"]

    z = [[matrix[g1][g2] for g2 in genes] for g1 in genes]
    fig = go.Figure(go.Heatmap(
        z=z, x=genes, y=genes,
        colorscale="YlOrRd",
        hovertemplate="<b>%{y}</b> + <b>%{x}</b><br>Co-occur: %{z} plasmids<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=550,
        xaxis=dict(tickangle=-45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=8), autorange="reversed"),
    )
    return fig


def make_integron_cassette_chart(integron_data):
    """Bar chart of gene cassettes found near qacEdelta1."""
    genes = integron_data.get("cassette_genes", [])
    if not genes:
        return go.Figure()

    # Normalize drug class names
    class_map = {
        "SULFONAMIDE": "Sulfonamide", "sulfonamide antibiotic": "Sulfonamide",
        "AMINOGLYCOSIDE": "Aminoglycoside", "aminoglycoside antibiotic": "Aminoglycoside",
        "TRIMETHOPRIM": "Trimethoprim", "diaminopyrimidine antibiotic": "Trimethoprim",
        "BETA-LACTAM": "Beta-lactam", "PHENICOL": "Phenicol",
        "RIFAMYCIN": "Rifamycin", "QUINOLONE": "Quinolone",
        "fluoroquinolone antibiotic": "Quinolone",
        "AMINOGLYCOSIDE/QUINOLONE": "Aminoglycoside/Quinolone",
        "fluoroquinolone antibiotic; aminoglycoside antibiotic": "Aminoglycoside/Quinolone",
        "QUATERNARY AMMONIUM": "QAC",
    }
    df = pd.DataFrame(genes)
    df["drug_class"] = df["drug_class"].fillna("").map(
        lambda x: class_map.get(x, x if x else "Other"))
    fig = px.bar(df, x="gene", y="cnt", color="drug_class",
                 color_discrete_sequence=px.colors.qualitative.Set2 + px.colors.qualitative.Pastel)
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis_tickangle=-45, xaxis_title="",
        yaxis_title="Plasmids with gene within 5kb of qacEdelta1",
        legend=dict(title="Drug Class", font_size=9),
    )
    return fig


def make_integron_mobility_chart(integron_data):
    """Donut: mobility of class 1 integron-carrying plasmids."""
    mob = integron_data.get("mobility", {})
    if not mob:
        return go.Figure()
    labels = [k.capitalize() for k in mob.keys()]
    values = list(mob.values())
    colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]
    fig = go.Figure(go.Pie(
        labels=labels, values=values, hole=0.55,
        marker_colors=colors[:len(labels)],
        textinfo="label+percent", textfont_size=12,
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=350, showlegend=False,
    )
    return fig


def make_comobilization_chart(comob_data):
    """Heatmap: conjugative Inc x mobilizable Inc co-occurrence."""
    pairs = comob_data.get("top_pairs", [])
    if not pairs:
        return go.Figure()

    # Build matrix
    conj_incs = sorted(set(c for c, m, n in pairs))[:10]
    mob_incs = sorted(set(m for c, m, n in pairs))[:10]
    pair_map = {(c, m): n for c, m, n in pairs}

    z = [[pair_map.get((c, m), 0) for m in mob_incs] for c in conj_incs]
    fig = go.Figure(go.Heatmap(
        z=z, x=[f"mob:{m}" for m in mob_incs],
        y=[f"conj:{c}" for c in conj_incs],
        colorscale="Purples", texttemplate="%{z}", textfont={"size": 9},
        hovertemplate="Conjugative: <b>%{y}</b><br>Mobilizable: <b>%{x}</b><br>Co-located: %{z}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=450,
        xaxis=dict(tickangle=-45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
    )
    return fig


# ---------------------------------------------------------------------------
def make_is_family_chart():
    """Bar chart of IS element family distribution."""
    families = db.is_family_counts()
    if not families:
        return go.Figure()
    # Top 15 families
    top = list(families.items())[:15]
    df = pd.DataFrame([
        {"Family": f, "Annotations": d["count"], "Plasmids": d["plasmids"]}
        for f, d in top
    ])
    fig = px.bar(df, x="Family", y="Annotations",
                 color_discrete_sequence=[COLORS["accent2"]],
                 hover_data={"Plasmids": True})
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=420,
        xaxis_tickangle=-45, xaxis_title="", yaxis_title="Total annotations",
    )
    return fig


def make_is_mobility_chart():
    """Stacked bar: IS families by mobility type."""
    fam_mob = db.is_family_by_mobility()
    if not fam_mob:
        return go.Figure()
    top_fams = sorted(fam_mob.keys(),
                      key=lambda f: sum(fam_mob[f].values()), reverse=True)[:12]
    mob_types = ["conjugative", "mobilizable", "non-mobilizable"]
    mob_colors = [COLORS["accent3"], COLORS["accent"], COLORS["accent4"]]

    fig = go.Figure()
    for mob, color in zip(mob_types, mob_colors):
        vals = [fam_mob[f].get(mob, 0) for f in top_fams]
        fig.add_trace(go.Bar(x=top_fams, y=vals, name=mob.capitalize(),
                             marker_color=color))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=420,
        barmode="stack", xaxis_tickangle=-45, xaxis_title="",
        yaxis_title="Annotations",
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1, font_size=11),
    )
    return fig


def make_is_inc_heatmap():
    """Heatmap: IS families vs Inc groups."""
    inc_is = db.is_family_by_inc()
    if not inc_is:
        return go.Figure()
    # Get top IS families and Inc groups
    all_is = {}
    for inc, fams in inc_is.items():
        for fam, cnt in fams.items():
            all_is[fam] = all_is.get(fam, 0) + cnt
    top_is = sorted(all_is, key=all_is.get, reverse=True)[:12]
    top_inc = sorted(inc_is.keys(),
                     key=lambda i: sum(inc_is[i].values()), reverse=True)[:10]

    z = [[inc_is[inc].get(fam, 0) for fam in top_is] for inc in top_inc]
    fig = go.Figure(go.Heatmap(
        z=z, x=top_is, y=top_inc,
        colorscale="Viridis", texttemplate="%{z}", textfont={"size": 9},
        hovertemplate="Inc: <b>%{y}</b><br>IS: <b>%{x}</b><br>Count: %{z:,}<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=420,
        xaxis=dict(tickangle=-45), yaxis=dict(autorange="reversed"),
    )
    return fig


def make_relaxase_t4ss_chart():
    """Heatmap: observed relaxase x T4SS co-occurrence with known rules overlay."""
    data, all_relax, all_mpf = db.relaxase_t4ss_compatibility()
    if not data:
        return go.Figure()

    # Build matrix
    matrix = {r: {m: 0 for m in all_mpf} for r in all_relax}
    known = {r: {m: False for m in all_mpf} for r in all_relax}
    for d in data:
        matrix[d["relaxase"]][d["mpf"]] = d["count"]
        known[d["relaxase"]][d["mpf"]] = d["known_compatible"]

    z = [[matrix[r][m] for m in all_mpf] for r in all_relax]
    text = []
    for r in all_relax:
        row = []
        for m in all_mpf:
            cnt = matrix[r][m]
            tag = " *" if known[r][m] else ""
            row.append(f"{cnt:,}{tag}")
        text.append(row)

    fig = go.Figure(go.Heatmap(
        z=z, x=[f"T4SS: {m}" for m in all_mpf],
        y=[f"Relaxase: {r}" for r in all_relax],
        colorscale="YlOrRd", text=text, texttemplate="%{text}",
        textfont={"size": 10},
        hovertemplate="<b>%{y}</b> + <b>%{x}</b><br>Co-located: %{z:,}<br>(* = known compatible)<extra></extra>",
    ))
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        yaxis=dict(autorange="reversed"),
    )
    return fig


def make_comobilization_ml_chart():
    """Random Forest to predict co-mobilization from plasmid features."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import LabelEncoder
    from sklearn.model_selection import cross_val_score
    import numpy as np

    positives, negatives = db.comobilization_ml_data()
    if len(positives) < 50 or len(negatives) < 50:
        return go.Figure(), None

    # Build feature matrix
    rows = []
    for r in positives:
        rows.append({**r, "colocated": 1})
    for r in negatives:
        r["conj_mpf"] = ""
        r["conj_inc"] = ""
        rows.append({**r, "colocated": 0})

    df = pd.DataFrame(rows)

    # Features
    features = []
    X_data = {}

    # Relaxase type (key feature)
    for relax in ["MOBF", "MOBP", "MOBQ", "MOBH", "MOBC", "MOBV"]:
        col = f"relax_{relax}"
        X_data[col] = df["relaxase_type"].fillna("").str.contains(relax).astype(int).values
        features.append(col)

    # OriT type
    for orit in ["MOBF", "MOBP", "MOBQ"]:
        col = f"orit_{orit}"
        X_data[col] = df["orit_type"].fillna("").str.contains(orit).astype(int).values
        features.append(col)

    # Plasmid size
    X_data["log_length"] = np.log1p(df["length"].fillna(0).values.astype(float))
    features.append("log_length")
    X_data["gc"] = df["gc"].fillna(0).values.astype(float)
    features.append("gc")

    # Genus
    le_genus = LabelEncoder()
    X_data["genus"] = le_genus.fit_transform(df["genus"].fillna("Unknown"))
    features.append("genus")

    # Inc groups
    for inc in ["IncFIB", "IncFII", "IncFIA", "IncI1", "IncN", "ColRNAI"]:
        col = f"inc_{inc}"
        X_data[col] = df["mob_inc"].fillna("").str.contains(inc).astype(int).values
        features.append(col)

    X = pd.DataFrame(X_data)[features]
    y = df["colocated"].values

    # Train with cross-validation
    rf = RandomForestClassifier(n_estimators=100, max_depth=8, random_state=42)
    cv_scores = cross_val_score(rf, X, y, cv=5, scoring="accuracy")
    rf.fit(X, y)

    # Feature importance
    importances = rf.feature_importances_
    display_names = [f.replace("relax_", "Relaxase: ").replace("orit_", "OriT: ")
                      .replace("inc_", "Inc: ").replace("log_", "Log ")
                      .replace("_", " ").capitalize()
                     for f in features]

    imp_df = pd.DataFrame({
        "Feature": display_names, "Importance": importances,
    }).sort_values("Importance", ascending=True)

    fig = px.bar(imp_df, x="Importance", y="Feature", orientation="h",
                 color="Importance", color_continuous_scale="Viridis")
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=500,
        coloraxis_showscale=False,
        xaxis_title="Feature Importance (Gini)", yaxis_title="",
    )

    accuracy = round(cv_scores.mean() * 100, 1)
    return fig, accuracy


# ---------------------------------------------------------------------------
# Charts will be pre-computed after all functions are defined (see below)
_shap_fig = _shap_acc = None

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

HEADER = html.Div(className="header", children=[
    html.Div(className="header-content", children=[
        html.Div(className="header-brand", children=[
            html.Span("Plasmid", className="brand-plasmid"),
            html.Span("NET", className="brand-net"),
        ]),
        html.Div(className="header-links", children=[
            html.A("Docs", href="https://plasmidnet.readthedocs.io/en/latest/",
                    target="_blank", className="header-link", title="Documentation"),
            html.A(
                html.Img(src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/github/github-original.svg",
                         style={"height": "20px", "width": "20px", "opacity": "0.7"}),
                href="https://github.com/lcerdeira/plasmidnet",
                target="_blank", className="header-link", title="View on GitHub",
            ),
        ]),
    ]),
])

_total_plasmids = overview.get("total", 0) + overview.get("ncbi_extra", 0)
STATS_ROW = html.Div(className="hero-banner", children=[
    html.Div(className="hero-inner", children=[
        html.Img(src="/assets/logo.png", className="hero-logo"),
        html.Div(className="hero-stats", children=[
            html.H2(f"{_total_plasmids:,}", className="hero-number"),
            html.P("Total Plasmids", className="hero-title"),
            html.Div(className="hero-breakdown", children=[
                html.Span(f"PLSDB {overview.get('total',0):,}", className="hero-tag"),
                html.Span("+", className="hero-plus"),
                html.Span(f"NCBI {overview.get('ncbi_extra',0):,}", className="hero-tag"),
            ]),
            html.Div(className="hero-badges", children=[
                html.Span(f"{overview.get('total_amr',0):,} AMR", className="badge"),
                html.Span(f"{overview.get('total_virulence_factors',0):,} Virulence",
                          className="badge"),
                html.Span(f"{overview.get('source_RefSeq',0):,} RefSeq", className="badge"),
                html.Span(f"{overview.get('source_INSDC',0):,} INSDC", className="badge"),
            ]),
        ]),
    ]),
    # About section
    html.Div(className="hero-about", children=[
        html.Details([
            html.Summary("About PlasmidNet", className="about-toggle"),
            html.Div(className="about-text", children=[
                html.P([
                    "Plasmids are the primary vehicles of antimicrobial resistance (AMR) "
                    "gene dissemination among bacteria, yet understanding their global "
                    "distribution, mobility, and evolutionary dynamics remains fragmented "
                    "across disconnected databases and tools."
                ]),
                html.P([
                    "PlasmidNet was created to bridge this gap \u2014 providing a single, "
                    "interactive platform that unifies ",
                    html.B("182,362 complete plasmids"),
                    " from PLSDB and NCBI with comprehensive analytics for AMR, "
                    "virulence, mobilization, geographic spread, and sequence analysis. "
                    "By combining classical epidemiological approaches with machine learning, "
                    "PlasmidNet enables researchers to ask questions that span the full "
                    "complexity of plasmid biology:"
                ]),
                html.Ul([
                    html.Li("Which incompatibility groups carry the most dangerous "
                            "resistance genes, and are they spreading?"),
                    html.Li("Are regional differences in plasmid mobility driven by "
                            "biology or sequencing bias?"),
                    html.Li("Which mobilizable plasmids are hitchhiking on conjugative "
                            "ones, and can we predict co-mobilization?"),
                    html.Li("Is a plasmid sequence natural or engineered \u2014 and what "
                            "cloning artifacts reveal about its history?"),
                ]),
                html.P([
                    "Built for the One Health community \u2014 from clinical microbiologists "
                    "tracking hospital outbreaks to environmental scientists monitoring "
                    "resistance in soil and water. PlasmidNet is open-source and freely "
                    "available, with a REST API for programmatic access and full "
                    "documentation at ",
                    html.A("plasmidnet.readthedocs.io",
                           href="https://plasmidnet.readthedocs.io/en/latest/",
                           target="_blank", style={"color": COLORS["accent"]}),
                    ".",
                ]),
            ]),
        ]),
    ]),
])

TABS = html.Div(className="tabs-container", children=[
    dcc.Tabs(
        id="main-tabs",
        value="overview",
        children=[
            dcc.Tab(label="Overview", value="overview", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Taxonomy", value="taxonomy", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="AMR", value="amr", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Viewer", value="viewer", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Inc/MOB/IS", value="mobility", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Correlations", value="correlations", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Geography", value="geography", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Analytics", value="analytics", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Compare", value="compare", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Seq Analysis", value="seqanalysis", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="API", value="api_explorer", className="tab", selected_className="tab--selected"),
        ],
    ),
])


def overview_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="chart-grid-3", children=[
            html.Div(className="chart-card", children=[
                html.H3("Topology Distribution", className="chart-title"),
                dcc.Graph(figure=make_topology_chart(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Data Source", className="chart-title"),
                dcc.Graph(figure=make_source_chart(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Host Kingdom", className="chart-title"),
                dcc.Graph(figure=make_kingdom_chart(), config={"displayModeBar": False}),
            ]),
        ]),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Temporal Distribution", className="chart-title"),
                dcc.Graph(figure=make_temporal_chart(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("GC Content Distribution", className="chart-title"),
                dcc.Graph(figure=make_gc_chart(), config={"displayModeBar": False}),
            ]),
        ]),
        html.Div(className="chart-grid-1", children=[
            html.Div(className="chart-card", children=[
                html.H3("Plasmid Length Distribution", className="chart-title"),
                dcc.Graph(figure=make_length_chart(), config={"displayModeBar": False}),
            ]),
        ]),
    ])


def taxonomy_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="chart-grid-1", children=[
            html.Div(className="chart-card", children=[
                html.H3("Top 20 Host Genera", className="chart-title"),
                html.P("Number of plasmids per bacterial genus",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_taxonomy_chart(), config={"displayModeBar": False}),
            ]),
        ]),
        html.Div(className="chart-card", style={"marginTop": "20px"}, children=[
            html.H3("Explore by Taxonomy", className="chart-title"),
            html.Div(className="filter-row", children=[
                html.Div([
                    html.Label("Genus", className="filter-label"),
                    dcc.Input(
                        id="taxonomy-genus-input", type="text",
                        placeholder="e.g., Escherichia", className="filter-input",
                    ),
                ]),
                html.Div([
                    html.Label("Species", className="filter-label"),
                    dcc.Input(
                        id="taxonomy-species-input", type="text",
                        placeholder="e.g., Escherichia coli", className="filter-input",
                    ),
                ]),
                html.Button("Search", id="taxonomy-search-btn", className="btn-primary"),
            ]),
            html.Div(id="taxonomy-search-result", className="search-result"),
        ]),
    ])


def amr_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("AMR Genes by Drug Class", className="chart-title"),
                html.P("Aggregated resistance gene counts by antimicrobial category",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_drug_class_summary(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Individual AMR Gene Frequency", className="chart-title"),
                html.P("Number of plasmids carrying each resistance gene",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_amr_chart(), config={"displayModeBar": False}),
            ]),
        ]),
    ])


def viewer_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="chart-card", children=[
            html.H3("Plasmid Architecture Viewer", className="chart-title"),
            html.P(
                "Enter an NCBI accession to visualize the plasmid architecture. "
                "Hover over features to see gene details.",
                className="chart-subtitle",
            ),
            html.Div(className="lookup-row", children=[
                dcc.Input(
                    id="viewer-input", type="text",
                    placeholder="e.g., NZ_CP031107.1",
                    className="filter-input lookup-input-wide",
                    value="",
                ),
                html.Button("Visualize", id="viewer-btn", className="btn-primary"),
            ]),
            html.P(
                "CDS genes, GC content, and AMR/MOB/BGC features are fetched "
                "from NCBI GenBank and PLSDB automatically.",
                className="chart-subtitle", style={"marginTop": "8px"},
            ),
            dcc.Loading(
                type="circle", color=COLORS["accent"],
                children=html.Div(id="viewer-result"),
            ),
        ]),
    ])


def mobility_tab():
    sample_note = ""
    ss = mob_typing.get("sample_size", 0)
    if ss and ss < 72000:
        sample_note = f" (based on sample of {ss:,} plasmids)"

    return html.Div(className="tab-content", children=[
        # Row 1: Inc groups + Mobility
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Incompatibility Groups (Replicon Types)", className="chart-title"),
                html.P(
                    f"Top 20 replicon type classifications{sample_note}",
                    className="chart-subtitle",
                ),
                dcc.Graph(figure=make_inc_group_chart(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Predicted Mobility", className="chart-title"),
                html.P(
                    "Distribution of conjugative, mobilizable, and non-mobilizable plasmids",
                    className="chart-subtitle",
                ),
                dcc.Graph(figure=make_mobility_chart(), config={"displayModeBar": False}),
            ]),
        ]),
        # Row 2: Relaxase + MPF
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Relaxase Types", className="chart-title"),
                html.P(
                    "Relaxase classification determines the mobilization mechanism",
                    className="chart-subtitle",
                ),
                dcc.Graph(figure=make_relaxase_chart(), config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("MPF Types (Mating Pair Formation)", className="chart-title"),
                html.P(
                    "Type IV secretion system classification for conjugative transfer",
                    className="chart-subtitle",
                ),
                dcc.Graph(figure=make_mpf_chart(), config={"displayModeBar": False}),
            ]),
        ]),

        # ── IS Elements & Transposons ──────────────────────────────
        html.Div(className="correlation-section", style={"marginTop": "28px"}, children=[
            html.H3("IS Elements & Transposon Families",
                     className="correlation-heading"),
            html.P("438,616 PGAP annotations across 72,556 plasmids",
                   className="chart-subtitle"),
        ]),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("IS Family Distribution", className="chart-title"),
                dcc.Graph(figure=make_is_family_chart(),
                          config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("IS Families by Mobility", className="chart-title"),
                dcc.Graph(figure=make_is_mobility_chart(),
                          config={"displayModeBar": False}),
            ]),
        ]),
        html.Div(className="chart-grid-1", children=[
            html.Div(className="chart-card", children=[
                html.H3("IS Families by Inc Group", className="chart-title"),
                dcc.Graph(figure=make_is_inc_heatmap(),
                          config={"displayModeBar": False}),
            ]),
        ]),
    ])


def correlations_tab():
    ss = correlations.get("sample_size", 0)
    note = f" (n = {ss:,} sampled plasmids)" if ss else ""

    return html.Div(className="tab-content", children=[
        # Section header
        html.Div(className="section-header", children=[
            html.H2("Cross-Feature Correlations", className="section-title"),
            html.P(
                f"Co-occurrence analysis across AMR, Inc groups, mobility, "
                f"heavy metals, phage elements, and pMLST{note}",
                className="chart-subtitle",
            ),
        ]),

        # --- 1. AMR vs Inc Groups & Mobility ---
        html.Div(className="correlation-section", children=[
            html.H3("1. AMR Drug Class vs Inc Groups & Mobility",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("AMR Drug Class x Inc Group Co-occurrence",
                             className="chart-title"),
                    html.P("Heatmap showing how resistance classes cluster with "
                           "specific incompatibility groups",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_amr_inc_heatmap(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("AMR Drug Class x Mobility",
                             className="chart-title"),
                    html.P("Resistance gene carriage by plasmid transfer mechanism",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_amr_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 2. Heavy Metal Resistance ---
        html.Div(className="correlation-section", children=[
            html.H3("2. Heavy Metal Resistance Genes",
                     className="correlation-heading"),
            html.Div(className="chart-grid-1", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Heavy Metal Gene Frequency", className="chart-title"),
                    html.P("Mercury, arsenic, copper, silver, and tellurium "
                           "resistance genes on plasmids",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_heavy_metal_genes_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Heavy Metal x Inc Group", className="chart-title"),
                    html.P("Which Inc groups carry metal resistance",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_heavy_metal_inc_heatmap(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Heavy Metal x Mobility", className="chart-title"),
                    html.P("Transfer potential of metal-resistant plasmids",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_heavy_metal_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 3. Phage & Mobile Genetic Elements ---
        html.Div(className="correlation-section", children=[
            html.H3("3. Phage & Mobile Genetic Elements",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Element Distribution", className="chart-title"),
                    html.P("Phage/prophage elements and transposases per plasmid",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_phage_distribution_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Phage Elements x Mobility", className="chart-title"),
                    html.P("Mean phage elements by mobility class",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_phage_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 4. pMLST Distribution ---
        html.Div(className="correlation-section", children=[
            html.H3("4. pMLST Distribution", className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("pMLST Schemes", className="chart-title"),
                    html.P("Plasmid MLST scheme frequency",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_pmlst_scheme_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("pMLST Allele Frequency", className="chart-title"),
                    html.P("Top 20 alleles across all typed plasmids",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_pmlst_alleles_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-1", children=[
                html.Div(className="chart-card", children=[
                    html.H3("pMLST Scheme x Mobility", className="chart-title"),
                    html.P("Transfer potential by plasmid MLST scheme",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_pmlst_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 5. Virulence Factors ---
        html.Div(className="correlation-section", children=[
            html.H3(
                f"5. Virulence Factors ({correlations.get('vir_total', 0):,} "
                f"plasmids with virulence genes)",
                className="correlation-heading",
            ),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Top Virulence Genes", className="chart-title"),
                    html.P("Most frequent virulence-associated genes",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_gene_frequency_chart(
                        "vir_gene_counts", "virulence", "#c0392b"),
                        config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Virulence x Mobility", className="chart-title"),
                    html.P("Transfer potential of virulence-carrying plasmids",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_mobility_chart(
                        "vir_vs_mobility", "virulence"),
                        config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-1", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Virulence x Inc Group", className="chart-title"),
                    html.P("Which Inc groups carry virulence factors",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_inc_chart(
                        "vir_vs_inc", "virulence", "#c0392b"),
                        config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 6. Toxin-Antitoxin Systems ---
        html.Div(className="correlation-section", children=[
            html.H3(
                f"6. Toxin-Antitoxin Systems ({correlations.get('ta_total', 0):,} "
                f"plasmids with TA genes)",
                className="correlation-heading",
            ),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Top TA Genes", className="chart-title"),
                    html.P("Most frequent toxin-antitoxin components (from PGAP)",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_gene_frequency_chart(
                        "ta_gene_counts", "TA", "#e67e22"),
                        config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("TA x Mobility", className="chart-title"),
                    html.P("Mobility of plasmids carrying TA systems",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_mobility_chart(
                        "ta_vs_mobility", "TA"),
                        config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-1", children=[
                html.Div(className="chart-card", children=[
                    html.H3("TA x Inc Group", className="chart-title"),
                    html.P("Which Inc groups carry TA systems",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_inc_chart(
                        "ta_vs_inc", "TA", "#e67e22"),
                        config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 7. Quaternary Ammonium Compounds ---
        html.Div(className="correlation-section", children=[
            html.H3(
                f"7. QAC Resistance ({correlations.get('qac_total', 0):,} "
                f"plasmids with QAC genes)",
                className="correlation-heading",
            ),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("QAC Resistance Genes", className="chart-title"),
                    html.P("Quaternary ammonium compound efflux genes",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_gene_frequency_chart(
                        "qac_gene_counts", "QAC", "#1abc9c"),
                        config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("QAC x Mobility", className="chart-title"),
                    html.P("Transfer potential of QAC-resistant plasmids",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_mobility_chart(
                        "qac_vs_mobility", "QAC"),
                        config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-1", children=[
                html.Div(className="chart-card", children=[
                    html.H3("QAC x Inc Group", className="chart-title"),
                    html.P("Which Inc groups carry QAC resistance",
                           className="chart-subtitle"),
                    dcc.Graph(figure=_make_feature_inc_chart(
                        "qac_vs_inc", "QAC", "#1abc9c"),
                        config={"displayModeBar": False}),
                ]),
            ]),
        ]),
    ])


def geography_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("Geographic Distribution", className="section-title"),
            html.P("57,751 geolocated plasmids from PLSDB across 100+ countries",
                   className="chart-subtitle"),
        ]),

        # --- 1. Global map ---
        html.Div(className="correlation-section", children=[
            html.H3("1. Global Plasmid Map by Mobility",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.P("Each dot is a plasmid coloured by predicted mobility. "
                       "Hover to see accession, country, and Inc group.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_global_map(),
                          config={"scrollZoom": True}),
            ]),
        ]),

        # --- 2. Country comparison ---
        html.Div(className="correlation-section", children=[
            html.H3("2. Country Comparison",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Mobility by Country (Top 20)",
                             className="chart-title"),
                    html.P("Conjugative, mobilizable, and non-mobilizable "
                           "plasmid counts per country",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_country_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Inc Groups by Country",
                             className="chart-title"),
                    html.P("Incompatibility group distribution across "
                           "top reporting countries",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_country_inc_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 3. Temporal animations ---
        html.Div(className="correlation-section", children=[
            html.H3("3. Temporal Spread Animations",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.H3("All Plasmids", className="chart-title"),
                html.P("Cumulative plasmid reports per country (2000-2024). "
                       "Press play to watch the global spread.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_temporal_animation(),
                          config={"scrollZoom": True}),
            ]),
            html.Div(className="chart-grid-2", style={"marginTop": "16px"}, children=[
                html.Div(className="chart-card", children=[
                    html.H3("Inc Group Spread", className="chart-title"),
                    html.P("Select an Inc group to see its temporal spread",
                           className="chart-subtitle"),
                    dcc.Dropdown(
                        id="geo-inc-dropdown",
                        options=[{"label": g, "value": g}
                                 for g in list(mob_typing["rep_types"].keys())[:15]],
                        value="IncFIB",
                        className="filter-input",
                        style={"width": "200px", "marginBottom": "10px"},
                    ),
                    dcc.Loading(html.Div(id="geo-inc-animation")),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Mobility Type Spread", className="chart-title"),
                    html.P("Select a mobility type to see its temporal spread",
                           className="chart-subtitle"),
                    dcc.Dropdown(
                        id="geo-mob-dropdown",
                        options=[{"label": m.capitalize(), "value": m}
                                 for m in ["conjugative", "mobilizable", "non-mobilizable"]],
                        value="conjugative",
                        className="filter-input",
                        style={"width": "200px", "marginBottom": "10px"},
                    ),
                    dcc.Loading(html.Div(id="geo-mob-animation")),
                ]),
            ]),
        ]),

        # --- 4. Host / Source correlations ---
        html.Div(className="correlation-section", children=[
            html.H3("4. Plasmid Host & Source",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Host / Source Distribution", className="chart-title"),
                    html.P("Classification: Human, Animal, Soil, Water, Food, Environment",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_host_distribution_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Mobility by Host", className="chart-title"),
                    html.P("Transfer potential varies significantly by source",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_host_mobility_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Inc Groups by Host", className="chart-title"),
                    html.P("Which replicon types dominate in each host category",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_host_inc_heatmap(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("AMR Drug Classes by Host", className="chart-title"),
                    html.P("Resistance profiles differ between human, animal, "
                           "and environmental sources",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_host_amr_heatmap(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),
    ])


def analytics_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("Statistical Analytics", className="section-title"),
            html.P("Confound decomposition: are regional mobility differences "
                   "biological or sampling bias?",
                   className="chart-subtitle"),
        ]),

        # --- 0. UMAP Clustering ---
        html.Div(className="correlation-section", children=[
            html.H3("Plasmid Feature Space (UMAP)", className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("Coloured by Mobility", className="chart-title"),
                    html.P("15,000 plasmids embedded in 2D by Inc group, genus, "
                           "host, length, GC. Clusters = natural plasmid ecotypes.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_umap_chart("mobility"),
                              config={"scrollZoom": True}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Coloured by Host Genus", className="chart-title"),
                    html.P("Same embedding coloured by bacterial genus. "
                           "Species-specific clusters reveal host-plasmid associations.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_umap_chart("genus"),
                              config={"scrollZoom": True}),
                ]),
            ]),
        ]),

        # --- 1. Matched comparison ---
        html.Div(className="correlation-section", children=[
            html.H3("1. Matched Comparison (Species-Controlled)",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.H3("Conjugative % by Species x Country", className="chart-title"),
                html.P("Each cell shows the conjugative fraction for a single species "
                       "in a single country. Green = more conjugative, red = less. "
                       "If regional differences persist within the same species, "
                       "the effect is not purely a species-composition artifact.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_matched_comparison(),
                          config={"displayModeBar": False}),
            ]),
        ]),

        # --- 2. Rarefaction ---
        html.Div(className="correlation-section", children=[
            html.H3("2. Rarefaction Analysis (Sample Size Effect)",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.H3("Conjugative Fraction vs Sample Size",
                         className="chart-title"),
                html.P("Each country is subsampled to equal n (50x bootstrap). "
                       "Error bars show standard deviation. If curves converge, "
                       "the differences are real. If they diverge, they may be noise. "
                       "Wide error bars at small n indicate insufficient data.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_rarefaction_chart(),
                          config={"displayModeBar": False}),
            ]),
        ]),

        # --- 3. XGBoost + SHAP Feature Importance ---
        html.Div(className="correlation-section", children=[
            html.H3("3. XGBoost + SHAP Feature Importance",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("SHAP Values (XGBoost)", className="chart-title"),
                    html.P("XGBoost trained on 57K plasmids. SHAP values show "
                           "each feature's contribution to the mobility prediction.",
                           className="chart-subtitle"),
                    html.P(f"Model accuracy: {_shap_acc}%",
                           className="chart-subtitle",
                           style={"fontWeight": "600", "color": COLORS["accent"]}),
                    dcc.Graph(figure=_shap_fig,
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Random Forest Importance", className="chart-title"),
                    html.P("Gini importance from a separate Random Forest model "
                           "for comparison with SHAP.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_feature_importance_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 4. Temporal Trends ---
        html.Div(className="correlation-section", children=[
            html.H3("4. Temporal Trends (2010-2024)",
                     className="correlation-heading"),
            html.Div(className="chart-grid-2", children=[
                html.Div(className="chart-card", children=[
                    html.H3("AMR Drug Class Prevalence Over Time",
                             className="chart-title"),
                    html.P("Percentage of plasmids carrying each resistance class "
                           "per year. Rising trends suggest increasing selective pressure.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_temporal_trends_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("Inc Group Prevalence Over Time",
                             className="chart-title"),
                    html.P("Percentage of plasmids with each Inc group per year. "
                           "Shifts indicate changing plasmid populations.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_inc_trends_chart(),
                              config={"displayModeBar": False}),
                ]),
            ]),
        ]),

        # --- 5. Simpson's paradox ---
        html.Div(className="correlation-section", children=[
            html.H3("5. Simpson's Paradox Detector",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.H3("Inc Groups Where Overall Trend Reverses by Species",
                         className="chart-title"),
                html.P("Cases where an Inc group appears mostly conjugative overall, "
                       "but is mostly non-mobilizable within individual species "
                       "(or vice versa). This reveals confounding by species composition.",
                       className="chart-subtitle"),
                make_simpson_paradox_table(),
            ]),
        ]),

        # --- 6. AMR Co-occurrence Network ---
        html.Div(className="correlation-section", children=[
            html.H3("6. AMR Gene Co-occurrence Network",
                     className="correlation-heading"),
            html.Div(className="chart-card", children=[
                html.H3("Which Resistance Genes Travel Together",
                         className="chart-title"),
                html.P("Heatmap showing AMR genes co-occurring on the same plasmid "
                       "(min 100 co-occurrences). Hot spots reveal "
                       "multi-drug resistance cassettes.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_amr_cooccurrence_chart(),
                          config={"displayModeBar": False}),
            ]),
        ]),

        # --- 7. Integron & Gene Cassette Analysis ---
        _build_integron_section(),

        # --- 8. Co-mobilization Analysis ---
        _build_comob_section(),

        # --- 9. Retro-mobilization & HGT Routes ---
        _build_retromob_section(),

        # --- 10. KPC Transposon Context ---
        _build_kpc_section(),

        # --- 11. Integron ML + Transposon-AMR ---
        _build_integron_ml_section(),
    ])


def is_families_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("IS Elements & Transposon Families", className="section-title"),
            html.P("Distribution of insertion sequences and transposons across "
                   "438,616 PGAP annotations from 72,556 PLSDB plasmids",
                   className="chart-subtitle"),
        ]),

        # Row 1: Family distribution + Mobility
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("IS Family Distribution", className="chart-title"),
                html.P("Top 15 IS element and transposon families by annotation count",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_is_family_chart(),
                          config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("IS Families by Mobility", className="chart-title"),
                html.P("Which IS families are enriched in conjugative vs "
                       "non-mobilizable plasmids",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_is_mobility_chart(),
                          config={"displayModeBar": False}),
            ]),
        ]),

        # Row 2: IS x Inc heatmap
        html.Div(className="chart-grid-1", children=[
            html.Div(className="chart-card", children=[
                html.H3("IS Families by Inc Group", className="chart-title"),
                html.P("Which IS element families dominate in each "
                       "incompatibility group. Different plasmid backbones "
                       "harbour distinct mobile element repertoires.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_is_inc_heatmap(),
                          config={"displayModeBar": False}),
            ]),
        ]),
    ])


def compare_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("Plasmid Comparison", className="section-title"),
            html.P("Compare up to 10 plasmids using BLASTn alignment "
                   "and pyCirclize visualization.",
                   className="chart-subtitle"),
        ]),
        html.Div(className="chart-card", children=[
            html.H3("Input Plasmids", className="chart-title"),
            html.P("Enter NCBI accessions (one per line) or upload FASTA files.",
                   className="chart-subtitle"),
            # Accessions input
            dcc.Textarea(
                id="compare-accessions",
                placeholder="MH595533\nMH595534\nNC_011102.1",
                style={"width": "100%", "height": "100px", "marginTop": "8px",
                       "fontFamily": "monospace", "fontSize": "0.9rem",
                       "borderRadius": "8px", "border": "1px solid #e2ddd5",
                       "padding": "10px", "background": "#faf7f2"},
            ),
            # File upload
            html.Div(style={"marginTop": "12px"}, children=[
                dcc.Upload(
                    id="compare-upload",
                    children=html.Div([
                        "Drag & drop or ",
                        html.A("click to upload", style={"color": COLORS["accent"],
                                                          "fontWeight": "600"}),
                        " FASTA / GenBank files",
                    ]),
                    className="upload-area",
                    multiple=True,
                ),
            ]),
            html.Div(style={"marginTop": "12px", "display": "flex", "gap": "10px"}, children=[
                html.Button("Compare", id="compare-btn", className="btn-primary"),
            ]),
        ]),
        dcc.Loading(
            type="circle", color=COLORS["accent"],
            children=html.Div(id="compare-results"),
        ),
    ])


def api_explorer_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("API Explorer", className="section-title"),
            html.P("Interactive access to the PlasmidNet REST API. "
                   "Full docs at /api/ endpoint.",
                   className="chart-subtitle"),
        ]),
        html.Div(className="chart-card", children=[
            html.H3("Query the Database", className="chart-title"),
            html.Div(className="filter-row", children=[
                html.Div([
                    html.Label("Endpoint", className="filter-label"),
                    dcc.Dropdown(
                        id="api-endpoint",
                        options=[
                            {"label": "Search plasmids", "value": "search"},
                            {"label": "AMR by gene", "value": "amr"},
                            {"label": "AMR by drug class", "value": "amr_class"},
                            {"label": "Typing by Inc group", "value": "typing"},
                            {"label": "Taxonomy by genus", "value": "taxonomy"},
                            {"label": "Plasmid details", "value": "plasmid"},
                            {"label": "Database stats", "value": "stats"},
                        ],
                        value="search",
                        style={"width": "250px"},
                    ),
                ]),
                html.Div([
                    html.Label("Query", className="filter-label"),
                    dcc.Input(
                        id="api-query", type="text",
                        placeholder="e.g., KPC, Klebsiella, IncFIB...",
                        className="filter-input",
                        style={"width": "300px"},
                    ),
                ]),
                html.Button("Run", id="api-run-btn", className="btn-primary"),
            ]),
            html.P(id="api-url-display", className="chart-subtitle",
                   style={"marginTop": "8px", "fontFamily": "monospace"}),
        ]),
        dcc.Loading(
            type="circle", color=COLORS["accent"],
            children=html.Div(id="api-results"),
        ),
    ])


def seqanalysis_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="section-header", children=[
            html.H2("Sequence Analysis", className="section-title"),
            html.P("Scan a DNA sequence for restriction sites, cryptic promoters, "
                   "RBS, codon bias, vector signatures, and engineering indicators.",
                   className="chart-subtitle"),
        ]),
        html.Div(className="chart-card", children=[
            html.H3("Input Sequence", className="chart-title"),
            html.P("Fetch from NCBI, upload a file, or paste/edit directly.",
                   className="chart-subtitle"),
            # Row 1: NCBI fetch
            html.Div(className="filter-row", children=[
                dcc.Input(
                    id="seq-accession", type="text",
                    placeholder="NCBI accession (e.g., MH595533)",
                    className="filter-input",
                    style={"width": "250px"},
                ),
                html.Button("Fetch & Analyze", id="seq-fetch-btn",
                            className="btn-primary"),
            ]),
            # Row 2: File upload
            dcc.Upload(
                id="seq-upload",
                children=html.Div([
                    "Drag & drop or ",
                    html.A("click to upload", style={"color": COLORS["accent"],
                                                      "fontWeight": "600"}),
                    " a FASTA or GenBank file",
                ]),
                className="upload-area",
                style={"marginTop": "12px"},
            ),
            # Row 3: Editable sequence textarea
            html.P("Or paste / edit sequence below:",
                   className="chart-subtitle", style={"marginTop": "12px"}),
            dcc.Textarea(
                id="seq-input",
                placeholder=">my_plasmid\nATGCATGCATGC...",
                style={"width": "100%", "height": "150px",
                       "fontFamily": "monospace", "fontSize": "0.85rem",
                       "borderRadius": "8px", "border": "1px solid #e2ddd5",
                       "padding": "10px", "background": "#faf7f2"},
            ),
            html.Button("Analyze Sequence", id="seq-paste-btn",
                        className="btn-secondary", style={"marginTop": "10px"}),
        ]),
        dcc.Loading(
            type="circle", color=COLORS["accent"],
            children=html.Div(id="seq-results"),
        ),
    ])


def lookup_tab():
    return html.Div(className="tab-content", children=[
        html.Div(className="chart-card", children=[
            html.H3("Plasmid Lookup", className="chart-title"),
            html.P("Enter an NCBI accession to view plasmid details",
                   className="chart-subtitle"),
            html.Div(className="lookup-row", children=[
                dcc.Input(
                    id="lookup-input", type="text",
                    placeholder="e.g., NZ_CP031107.1",
                    className="filter-input lookup-input-wide",
                ),
                html.Button("Lookup", id="lookup-btn", className="btn-primary"),
            ]),
            dcc.Loading(
                id="lookup-loading",
                type="circle",
                color=COLORS["accent"],
                children=html.Div(id="lookup-result", className="lookup-result"),
            ),
        ]),
    ])


FOOTER = html.Div(className="footer", children=[
    html.P([
        "Copyright © 2026 PlasmidNet. All rights reserved. | ",
        html.A("", href="https://github.com/lcerdeira/plasmidnet",
               target="_blank"),
        " | Built with ",
        html.A("Dash", href="https://dash.plotly.com/", target="_blank"),
        " & ",
        html.A("Plotly", href="https://plotly.com/", target="_blank"),
    ]),
])

DOWNLOAD_BAR = html.Div(className="download-bar", children=[
    html.Div(className="download-bar-inner", children=[
        html.Span("Export Data:", className="download-label"),
        dcc.Dropdown(
            id="download-dataset",
            options=[{"label": label, "value": key}
                     for key, (label, _) in db.EXPORT_QUERIES.items()],
            value=None,
            placeholder="Select dataset...",
            className="download-dropdown",
            style={"width": "250px"},
        ),
        html.Button("Download CSV", id="download-btn", className="btn-secondary",
                     style={"marginLeft": "8px"}),
        dcc.Download(id="download-csv"),
    ]),
])

app.layout = html.Div(className="app-container", children=[
    HEADER,
    html.Div(className="main-content", children=[
        STATS_ROW,
        TABS,
        html.Div(id="tab-content"),
    ]),
    DOWNLOAD_BAR,
    FOOTER,
])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(Output("tab-content", "children"), Input("main-tabs", "value"))
def render_tab(tab):
    if tab == "overview":
        return overview_tab()
    elif tab == "taxonomy":
        return taxonomy_tab()
    elif tab == "amr":
        return amr_tab()
    elif tab == "viewer":
        return viewer_tab()
    elif tab == "mobility":
        return mobility_tab()
    elif tab == "correlations":
        return correlations_tab()
    elif tab == "geography":
        return geography_tab()
    elif tab == "analytics":
        return analytics_tab()
    elif tab == "compare":
        return compare_tab()
    elif tab == "api_explorer":
        return api_explorer_tab()
    elif tab == "seqanalysis":
        return seqanalysis_tab()
    elif tab == "lookup":
        return lookup_tab()
    return overview_tab()


@callback(
    Output("compare-results", "children"),
    Input("compare-btn", "n_clicks"),
    State("compare-accessions", "value"),
    State("compare-upload", "contents"),
    State("compare-upload", "filename"),
    prevent_initial_call=True,
)
def run_comparison(n_clicks, accessions_text, upload_contents, upload_filenames):
    import base64 as b64

    accessions = []
    if accessions_text:
        accessions = [a.strip() for a in accessions_text.strip().split("\n") if a.strip()]

    # Parse uploaded files
    uploaded_seqs = []  # list of (name, sequence)
    if upload_contents:
        if not isinstance(upload_contents, list):
            upload_contents = [upload_contents]
            upload_filenames = [upload_filenames]
        for content, fname in zip(upload_contents, upload_filenames):
            _, content_string = content.split(",")
            decoded = b64.b64decode(content_string).decode("utf-8", errors="ignore")
            # Parse FASTA
            seqs = _parse_multi_fasta(decoded)
            for name, seq in seqs:
                uploaded_seqs.append((fname.rsplit(".", 1)[0] if not name else name, seq))

    total = len(accessions) + len(uploaded_seqs)
    if total < 2:
        return html.P("Provide at least 2 plasmids (accessions + uploads).",
                      className="text-muted")
    if total > 10:
        return html.P("Maximum 10 plasmids.", className="text-muted")

    # Fetch NCBI accessions
    sequences = []
    names = []
    for acc in accessions[:10]:
        gb_text = _fetch_genbank_text(acc)
        if not gb_text:
            return html.P(f"Could not fetch {acc} from NCBI.", className="text-danger")
        seq = _parse_genbank_sequence(gb_text)
        if not seq or len(seq) < 100:
            return html.P(f"{acc}: no sequence data or too short.", className="text-danger")
        sequences.append(seq)
        names.append(acc)

    # Add uploaded sequences
    for name, seq in uploaded_seqs:
        if len(seq) >= 100:
            sequences.append(seq)
            names.append(name)

    if len(sequences) < 2:
        return html.P("Need at least 2 valid sequences.", className="text-muted")

    # Fetch CDS features for each sequence
    from data_loader import fetch_genbank_features
    all_features = []
    for name in names:
        feats = fetch_genbank_features(name)
        all_features.append(feats)

    # Run comparison with features
    img_b64, alignments = plasmid_compare.compare_plasmids(
        names, sequences, names, all_features)

    # Build results
    results = []

    # Summary badges
    badges = []
    for n, s, f in zip(names, sequences, all_features):
        label = f"{n}: {len(s):,} bp"
        if f:
            label += f" ({len(f)} CDS)"
        else:
            label += " (no annotations)"
        badges.append(html.Span(label, className="badge"))
    results.append(html.Div(className="feature-badges", children=badges))

    # Comparison image
    if img_b64:
        results.append(html.Img(src=img_b64, className="plasmid-map-img"))

    # Alignment table
    if alignments:
        results.append(html.Div(className="chart-card", style={"marginTop": "16px"}, children=[
            html.H3(f"BLASTn Alignments ({len(alignments)} HSPs)", className="chart-title"),
            html.Table(className="analytics-table", children=[
                html.Thead(html.Tr([
                    html.Th("Query"), html.Th("Subject"),
                    html.Th("Identity"), html.Th("Length"),
                    html.Th("Q start-end"), html.Th("S start-end"),
                ])),
                html.Tbody([
                    html.Tr([
                        html.Td(a["query"]), html.Td(a["subject"]),
                        html.Td(f"{a['identity']:.1f}%"),
                        html.Td(f"{a['length']:,}"),
                        html.Td(f"{a['q_start']:,}-{a['q_end']:,}"),
                        html.Td(f"{a['s_start']:,}-{a['s_end']:,}"),
                    ]) for a in alignments[:20]
                ]),
            ]),
        ]))
    elif not alignments:
        results.append(html.P("No BLASTn alignments found (sequences may be too divergent "
                              "or BLAST not installed on this server).",
                              className="text-muted", style={"marginTop": "12px"}))

    # CDS annotation tables for each plasmid
    for n, feats in zip(names, all_features):
        if feats:
            results.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
                html.H3(f"CDS Annotations: {n} ({len(feats)} genes)", className="chart-title"),
                html.Table(className="analytics-table", children=[
                    html.Thead(html.Tr([
                        html.Th("Gene"), html.Th("Product"),
                        html.Th("Start"), html.Th("End"), html.Th("Strand"),
                    ])),
                    html.Tbody([
                        html.Tr([
                            html.Td(f.get("gene") or f.get("locus_tag", ""),
                                    style={"fontWeight": "600"}),
                            html.Td(f.get("product", "")[:60]),
                            html.Td(f"{f['start']:,}"),
                            html.Td(f"{f['end']:,}"),
                            html.Td("+" if f.get("strand", 1) == 1 else "-"),
                        ]) for f in feats
                    ]),
                ]),
            ]))

    # Download section
    download_items = []
    for i, (n, s, f) in enumerate(zip(names, sequences, all_features)):
        fasta = plasmid_compare.export_fasta(n, s)
        gbk = plasmid_compare.export_genbank(n, s, f)
        download_items.append(html.Div(style={"marginTop": "8px"}, children=[
            html.Span(f"{n}: ", style={"fontWeight": "600", "marginRight": "8px"}),
            html.A("FASTA", href=f"data:text/plain,{fasta}",
                   download=f"{n}.fasta", className="btn-secondary",
                   style={"fontSize": "0.8rem", "marginRight": "6px", "padding": "4px 12px"}),
            html.A("GenBank", href=f"data:text/plain,{gbk}",
                   download=f"{n}.gbk", className="btn-secondary",
                   style={"fontSize": "0.8rem", "padding": "4px 12px"}),
        ]))

    results.append(html.Div(className="chart-card", style={"marginTop": "16px"}, children=[
        html.H3("Download Annotations", className="chart-title"),
        html.P("Download each plasmid as FASTA or annotated GenBank format.",
               className="chart-subtitle"),
    ] + download_items))

    return html.Div(results)


@callback(
    Output("api-results", "children"),
    Input("api-run-btn", "n_clicks"),
    State("api-endpoint", "value"),
    State("api-query", "value"),
    prevent_initial_call=True,
)
def run_api_query(n_clicks, endpoint, query):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate

    query = (query or "").strip()
    endpoint_map = {
        "search": ("search", {"q": query}),
        "amr": ("amr", {"gene": query}),
        "amr_class": ("amr", {"class": query}),
        "typing": ("typing", {"inc": query}),
        "taxonomy": ("taxonomy", {"genus": query}),
        "plasmid": (f"plasmid/{query}", {}),
        "stats": ("stats", {}),
    }

    path, params = endpoint_map.get(endpoint, ("", {}))
    if not path:
        return html.P("Invalid endpoint.", className="text-danger")

    # Query the API internally
    import json as json_mod
    if endpoint == "plasmid":
        data = db.plasmid_summary(query)
        if not data:
            return html.P(f"Accession '{query}' not found.", className="text-danger")
    elif endpoint == "stats":
        data = db.overview_stats()
        data["mobility"] = db.mobility_distribution()
    elif endpoint == "search":
        rows = db.q("""
            SELECT n.NUCCORE_ACC, n.NUCCORE_Description, n.NUCCORE_Length,
                   t.TAXONOMY_genus
            FROM nuccore n
            LEFT JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
            WHERE n.NUCCORE_ACC LIKE ? OR n.NUCCORE_Description LIKE ?
               OR t.TAXONOMY_genus LIKE ?
            LIMIT 50
        """, (f"%{query}%", f"%{query}%", f"%{query}%"))
        data = {"count": len(rows), "results": rows}
    elif endpoint == "amr":
        rows = db.q("SELECT NUCCORE_ACC, gene_symbol, drug_class FROM amr WHERE gene_symbol LIKE ? LIMIT 50",
                    (f"%{query}%",))
        data = {"count": len(rows), "results": rows}
    elif endpoint == "amr_class":
        rows = db.q("SELECT NUCCORE_ACC, gene_symbol, drug_class FROM amr WHERE UPPER(drug_class) LIKE ? LIMIT 50",
                    (f"%{query.upper()}%",))
        data = {"count": len(rows), "results": rows}
    elif endpoint == "typing":
        rows = db.q("SELECT NUCCORE_ACC, rep_type, predicted_mobility FROM typing WHERE rep_type LIKE ? LIMIT 50",
                    (f"%{query}%",))
        data = {"count": len(rows), "results": rows}
    elif endpoint == "taxonomy":
        rows = db.q("""
            SELECT n.NUCCORE_ACC, n.NUCCORE_Description, t.TAXONOMY_genus
            FROM nuccore n JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
            WHERE t.TAXONOMY_genus LIKE ? LIMIT 50
        """, (f"%{query}%",))
        data = {"count": len(rows), "results": rows}
    else:
        data = {}

    # Format as JSON + table
    json_str = json_mod.dumps(data, indent=2, default=str)

    results_children = [
        html.Pre(json_str[:3000], style={
            "background": "#faf7f2", "border": "1px solid #e2ddd5",
            "borderRadius": "8px", "padding": "12px", "fontSize": "0.8rem",
            "maxHeight": "400px", "overflow": "auto", "whiteSpace": "pre-wrap",
        }),
    ]

    # If results have a table, show it
    if isinstance(data, dict) and "results" in data and data["results"]:
        rows = data["results"][:20]
        if isinstance(rows[0], dict):
            headers = list(rows[0].keys())
            results_children.append(
                html.Table(className="analytics-table", style={"marginTop": "12px"}, children=[
                    html.Thead(html.Tr([html.Th(h) for h in headers])),
                    html.Tbody([
                        html.Tr([html.Td(str(r.get(h, ""))[:50]) for h in headers])
                        for r in rows
                    ]),
                ])
            )

    return html.Div(results_children)


@callback(
    Output("seq-input", "value"),
    Input("seq-upload", "contents"),
    State("seq-upload", "filename"),
    prevent_initial_call=True,
)
def load_uploaded_seq(contents, filename):
    """Load uploaded file content into the editable textarea."""
    if not contents:
        raise dash.exceptions.PreventUpdate
    import base64 as b64
    content_type, content_string = contents.split(",")
    decoded = b64.b64decode(content_string).decode("utf-8", errors="ignore")
    return decoded


@callback(
    Output("seq-results", "children"),
    Input("seq-fetch-btn", "n_clicks"),
    Input("seq-paste-btn", "n_clicks"),
    State("seq-accession", "value"),
    State("seq-input", "value"),
    prevent_initial_call=True,
)
def run_seq_analysis(fetch_clicks, paste_clicks, accession, pasted_seq):
    ctx = dash.callback_context
    triggered = ctx.triggered[0]["prop_id"] if ctx.triggered else ""

    seq = ""
    name = "query"

    if "fetch" in triggered and accession:
        accession = accession.strip()
        gb_text = _fetch_genbank_text(accession)
        if not gb_text:
            return html.P(f"Could not fetch {accession} from NCBI.",
                          className="text-danger")
        seq = _parse_genbank_sequence(gb_text)
        name = accession
    elif "paste" in triggered and pasted_seq:
        # Handle FASTA or raw
        text = pasted_seq.strip()
        if text.startswith(">"):
            lines = text.split("\n")
            name = lines[0][1:].strip() or "pasted"
            seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
        else:
            seq = text
        seq = "".join(c for c in seq.upper() if c in "ACGTN")
    else:
        return html.P("Enter an accession or paste a sequence.", className="text-muted")

    if len(seq) < 100:
        return html.P("Sequence too short (minimum 100 bp).", className="text-danger")

    # Run analysis
    results = seq_analysis.analyze_sequence(seq, name)
    if "error" in results:
        return html.P(results["error"], className="text-danger")

    # Build results display
    eng_score = results["engineering_score"]
    eng_color = "#10b981" if eng_score < 25 else "#f59e0b" if eng_score < 50 else "#ef4444"
    eng_label = results.get("classification", "Unknown")

    # Score card
    score_card = html.Div(className="chart-card", style={"marginTop": "16px"}, children=[
        html.Div(style={"display": "flex", "alignItems": "center", "gap": "20px"}, children=[
            html.Div(style={"textAlign": "center"}, children=[
                html.H2(f"{eng_score}", style={
                    "fontSize": "3rem", "fontWeight": "800", "color": eng_color,
                    "margin": "0", "lineHeight": "1"}),
                html.P("/100", style={"color": COLORS["text_muted"], "margin": "0"}),
                html.P(eng_label, style={
                    "fontWeight": "600", "color": eng_color, "margin": "4px 0 0 0"}),
            ]),
            html.Div([
                html.H3(f"{name} ({results['length']:,} bp, GC {results['gc_content']}%)",
                         className="chart-title"),
                html.Div(className="feature-badges", children=[
                    html.Span(f"{results['restriction_density']['total_sites']} restriction sites "
                              f"({results['restriction_density']['density_per_kb']}/kb)", className="badge"),
                    html.Span(f"{len(results['restriction_density']['hotspots'])} RE hotspots", className="badge"),
                    html.Span(f"{len(results['promoters'])} promoters", className="badge"),
                    html.Span(f"{len(results['rbs_sites'])} RBS sites", className="badge"),
                    html.Span(f"{len(results['vector_signatures'])} vector signatures", className="badge"),
                    html.Span(f"{len(results.get('is_elements',[]))} IS/Tn elements", className="badge"),
                    html.Span(f"{len(results.get('engineering_scars',[]))} eng. scars", className="badge"),
                ]),
            ]),
        ]),
    ])

    # Detail sections
    sections = []

    # Restriction sites
    re_sites = results["restriction_sites"]
    if re_sites:
        re_df = pd.DataFrame(re_sites[:30])
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"Restriction Sites ({len(re_sites)} total)", className="chart-title"),
            html.Table(className="analytics-table", children=[
                html.Thead(html.Tr([html.Th("Enzyme"), html.Th("Position"), html.Th("Sequence")])),
                html.Tbody([
                    html.Tr([html.Td(s["enzyme"]), html.Td(f"{s['position']:,}"),
                             html.Td(s["sequence"], style={"fontFamily": "monospace"})])
                    for s in re_sites[:20]
                ]),
            ]),
        ]))

    # Vector signatures
    if results["vector_signatures"]:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3("Vector Backbone Signatures Detected", className="chart-title"),
            html.P("These sequences match known cloning vector components.",
                   className="chart-subtitle", style={"color": "#ef4444"}),
            html.Ul([
                html.Li(f"{v['name']} at position {v['position']:,}")
                for v in results["vector_signatures"]
            ]),
        ]))

    # Promoters
    if results["promoters"]:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"Predicted Promoters ({len(results['promoters'])})", className="chart-title"),
            html.Table(className="analytics-table", children=[
                html.Thead(html.Tr([html.Th("Position"), html.Th("Strand"),
                                     html.Th("-35"), html.Th("-10"),
                                     html.Th("Spacer"), html.Th("Score")])),
                html.Tbody([
                    html.Tr([
                        html.Td(f"{p['position']:,}"), html.Td(p["strand"]),
                        html.Td(p["minus35"], style={"fontFamily": "monospace"}),
                        html.Td(p["minus10"], style={"fontFamily": "monospace"}),
                        html.Td(f"{p['spacer']} bp"), html.Td(f"{p['score']}/100"),
                    ]) for p in results["promoters"][:10]
                ]),
            ]),
        ]))

    # K-mer naturalness
    kmer = results["kmer_score"]
    sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
        html.H3("K-mer Naturalness Analysis", className="chart-title"),
        html.Div(className="detail-grid", children=[
            html.Div(className="detail-item", children=[
                html.Span("4-mer Entropy Ratio", className="detail-label"),
                html.Span(f"{kmer.get('entropy_ratio', 'N/A')}", className="detail-value"),
            ]),
            html.Div(className="detail-item", children=[
                html.Span("Palindrome Fraction", className="detail-label"),
                html.Span(f"{kmer.get('palindrome_fraction', 'N/A')}", className="detail-value"),
            ]),
            html.Div(className="detail-item", children=[
                html.Span("K-mer Score", className="detail-label"),
                html.Span(f"{kmer.get('score', 'N/A')}/100 ({kmer.get('note', '')})",
                          className="detail-value"),
            ]),
        ]),
    ]))

    # IS elements
    is_elements = results.get("is_elements", [])
    if is_elements:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"IS Elements & Transposons ({len(is_elements)} found)",
                     className="chart-title"),
            html.P("Natural mobile element markers. Their presence near RE hotspots "
                   "indicates natural transposon boundaries, not cloning junctions.",
                   className="chart-subtitle",
                   style={"color": COLORS["accent3"]}),
            html.Ul([html.Li(f"{e['name']} at position {e['position']:,}")
                     for e in is_elements[:10]]),
        ]))

    # Engineering scars
    eng_scars = results.get("engineering_scars", [])
    if eng_scars:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"Engineering Junction Scars ({len(eng_scars)} found)",
                     className="chart-title"),
            html.P("These patterns are specific to engineered constructs "
                   "(Golden Gate overhangs, Gibson Assembly junctions).",
                   className="chart-subtitle",
                   style={"color": COLORS["danger"]}),
            html.Ul([html.Li(f"{s['name']} at position {s['position']:,}: "
                             f"{s['sequence']}")
                     for s in eng_scars[:10]]),
        ]))

    # Score reasoning
    reasons = results.get("score_reasons", [])
    if reasons:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3("Score Reasoning", className="chart-title"),
            html.Ul([html.Li(r) for r in reasons]),
        ]))

    # NTEKPC detection
    ntekpc = results.get("ntekpc", [])
    if ntekpc:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"KPC Transposon Detection ({len(ntekpc)} found)",
                     className="chart-title"),
            html.Ul([
                html.Li([
                    html.Span(n["type"], style={"fontWeight": "700",
                              "color": COLORS["danger"]}),
                    f" at position {n['kpc_position']:,} — {n['note']}",
                ]) for n in ntekpc
            ]),
        ]))

    # Direct repeats (TSDs)
    dr = results.get("direct_repeats", [])
    if dr:
        sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
            html.H3(f"Direct Repeats / Target Site Duplications ({len(dr)} found)",
                     className="chart-title"),
            html.P("Short direct repeats flanking inserted elements — "
                   "hallmark of transposon insertion.",
                   className="chart-subtitle"),
            html.Table(className="analytics-table", children=[
                html.Thead(html.Tr([
                    html.Th("Sequence"), html.Th("Length"),
                    html.Th("Position 1"), html.Th("Position 2"),
                    html.Th("Distance"),
                ])),
                html.Tbody([
                    html.Tr([
                        html.Td(r["sequence"], style={"fontFamily": "monospace"}),
                        html.Td(f"{r['length']} bp"),
                        html.Td(f"{r['pos1']:,}"),
                        html.Td(f"{r['pos2']:,}"),
                        html.Td(f"{r['distance']:,} bp"),
                    ]) for r in dr[:10]
                ]),
            ]),
        ]))

    # Retro-mobilization (needs GenBank features)
    if "fetch" in triggered and accession:
        from data_loader import fetch_genbank_features
        feats = fetch_genbank_features(accession)
        if feats:
            retro = seq_analysis.detect_retromobilization(feats)
            cat_color = {
                "Self-transmissible (conjugative)": COLORS["accent3"],
                "Mobilizable": COLORS["accent"],
                "Retro-mobilizable (oriT only)": COLORS["accent4"],
                "Mobilizable (MOB genes)": COLORS["accent"],
                "Non-mobilizable": COLORS["text_muted"],
            }
            sections.append(html.Div(className="chart-card", style={"marginTop": "12px"}, children=[
                html.H3("Mobilization Assessment", className="chart-title"),
                html.P([
                    html.Span(retro["category"],
                              style={"fontWeight": "700", "fontSize": "1.1rem",
                                     "color": cat_color.get(retro["category"], COLORS["text"])}),
                ]),
                html.P(retro["detail"], className="chart-subtitle"),
                html.Div(className="detail-grid", children=[
                    html.Div(className="detail-item", children=[
                        html.Span("oriT", className="detail-label"),
                        html.Span("Yes" if retro["has_orit"] else "No", className="detail-value"),
                    ]),
                    html.Div(className="detail-item", children=[
                        html.Span("Relaxase", className="detail-label"),
                        html.Span("Yes" if retro["has_relaxase"] else "No", className="detail-value"),
                    ]),
                    html.Div(className="detail-item", children=[
                        html.Span("T4SS", className="detail-label"),
                        html.Span("Yes" if retro["has_t4ss"] else "No", className="detail-value"),
                    ]),
                    html.Div(className="detail-item", children=[
                        html.Span("MOB genes", className="detail-label"),
                        html.Span(", ".join(retro["mob_genes"][:5]) or "None",
                                  className="detail-value"),
                    ]),
                ]),
            ]))

    return html.Div([score_card] + sections)


@callback(
    Output("taxonomy-search-result", "children"),
    Input("taxonomy-search-btn", "n_clicks"),
    State("taxonomy-genus-input", "value"),
    State("taxonomy-species-input", "value"),
    prevent_initial_call=True,
)
def search_taxonomy(n_clicks, genus, species):
    if not genus and not species:
        return html.P("Please enter a genus or species name.", className="text-muted")

    try:
        where = []
        params = []
        if genus:
            where.append("t.TAXONOMY_genus LIKE ?")
            params.append(f"%{genus}%")
        if species:
            where.append("t.TAXONOMY_species LIKE ?")
            params.append(f"%{species}%")

        count = db.scalar(
            f"SELECT COUNT(*) FROM nuccore n JOIN taxonomy t "
            f"ON n.TAXONOMY_UID = t.TAXONOMY_UID WHERE {' AND '.join(where)}",
            tuple(params),
        )
        search_term = species or genus
        if count:
            return html.Div([
                html.P(f"Found {count:,} plasmids for {search_term}",
                       className="result-count"),
                html.A(
                    f"View in PLSDB",
                    href=f"https://ccb-microbe.cs.uni-saarland.de/plsdb2025/browse",
                    target="_blank", className="btn-secondary",
                ),
            ])
        else:
            return html.P(f"No plasmids found for {search_term}.", className="text-muted")
    except Exception as e:
        return html.P(f"Error: {str(e)}", className="text-danger")


@callback(
    Output("download-csv", "data"),
    Input("download-btn", "n_clicks"),
    State("download-dataset", "value"),
    prevent_initial_call=True,
)
def download_data(n_clicks, dataset):
    if not n_clicks or not dataset:
        raise dash.exceptions.PreventUpdate
    csv_str = db.export_csv(dataset)
    if not csv_str:
        raise dash.exceptions.PreventUpdate
    return dcc.send_string(csv_str, filename=f"plasmidnet_{dataset}.csv")




def _build_integron_section():
    """Build integron & gene cassette section."""
    integ = db.integron_analysis()
    return html.Div(className="correlation-section", children=[
        html.H3(f"7. Class 1 Integrons & Gene Cassettes "
                 f"({integ['class1_total']:,} plasmids)",
                 className="correlation-heading"),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Gene Cassettes near qacEdelta1", className="chart-title"),
                html.P("AMR genes within 5 kb of qacEdelta1 (3' conserved segment "
                       "of class 1 integrons). sul1 + aadA + dfrA is the classic pattern.",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_integron_cassette_chart(integ),
                          config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Integron Mobility", className="chart-title"),
                html.P(f"Mobility of {integ['class1_total']:,} class 1 "
                       f"integron-carrying plasmids (qacEdelta1 + sul1)",
                       className="chart-subtitle"),
                dcc.Graph(figure=make_integron_mobility_chart(integ),
                          config={"displayModeBar": False}),
            ]),
        ]),
    ])


def _build_comob_section():
    """Build co-mobilization section with relaxase-T4SS and ML."""
    comob = db.comobilization_data()
    ml_fig, ml_acc = make_comobilization_ml_chart()

    ml_section = []
    if ml_acc:
        ml_section = [
            html.Div(className="chart-grid-2", style={"marginTop": "16px"}, children=[
                html.Div(className="chart-card", children=[
                    html.H3("Relaxase x T4SS Compatibility", className="chart-title"),
                    html.P("Observed co-location of mobilizable relaxase types with "
                           "conjugative T4SS types. Values marked * are known "
                           "compatible pairs from the literature.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=make_relaxase_t4ss_chart(),
                              config={"displayModeBar": False}),
                ]),
                html.Div(className="chart-card", children=[
                    html.H3("ML: What Predicts Co-mobilization?",
                             className="chart-title"),
                    html.P(f"Random Forest (5-fold CV accuracy: {ml_acc}%) "
                           f"trained to predict whether a mobilizable plasmid "
                           f"shares a host with a conjugative one.",
                           className="chart-subtitle"),
                    dcc.Graph(figure=ml_fig, config={"displayModeBar": False}),
                ]),
            ]),
        ]

    return html.Div(className="correlation-section", children=[
        html.H3(f"8. Co-mobilization Analysis "
                 f"({comob['pct_colocated']}% of mobilizable plasmids "
                 f"share a host with a conjugative plasmid)",
                 className="correlation-heading"),
        html.Div(className="chart-card", children=[
            html.H3("Conjugative x Mobilizable Inc Group Pairs",
                     className="chart-title"),
            html.P("Which conjugative Inc groups co-exist with which mobilizable "
                   "Inc groups in the same bacterial host.",
                   className="chart-subtitle"),
            dcc.Graph(figure=make_comobilization_chart(comob),
                      config={"displayModeBar": False}),
        ]),
    ] + ml_section)


def _build_retromob_section():
    """Build retro-mobilization section (called directly, not via callback)."""
    retro = db.retromobilization_analysis()

    # Transfer routes bar chart
    routes = retro["routes"]
    route_df = pd.DataFrame([
        {"Route": k, "Plasmids": v} for k, v in routes.items()
    ]).sort_values("Plasmids", ascending=True)
    route_fig = px.bar(route_df, x="Plasmids", y="Route", orientation="h",
                       color_discrete_sequence=[COLORS["accent2"]])
    route_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=300,
        xaxis_title="Non-mobilizable plasmids", yaxis_title="",
    )

    # Retro-mob Inc groups bar
    inc_df = pd.DataFrame(
        sorted(retro["retro_inc"].items(), key=lambda x: -x[1])[:10],
        columns=["Inc Group (conjugative)", "Non-mob partners"],
    )
    inc_fig = px.bar(inc_df, x="Non-mob partners",
                     y="Inc Group (conjugative)", orientation="h",
                     color_discrete_sequence=[COLORS["accent4"]])
    inc_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=350,
        yaxis=dict(autorange="reversed"),
    )

    pct_retro = retro["retro_mob"] / retro["total_nonmob"] * 100

    return html.Div(className="correlation-section", children=[
        html.H3(f"9. Retro-mobilization & HGT Routes for Non-mobilizable Plasmids",
                 className="correlation-heading"),
        html.P([
            f"{retro['total_nonmob']:,} non-mobilizable plasmids lack autonomous "
            f"transfer machinery, yet ",
            html.B(f"{retro['retro_mob']:,} ({pct_retro:.0f}%)"),
            " share a bacterial host with a conjugative plasmid — enabling "
            "retro-mobilization (reverse transfer), where the conjugative "
            "plasmid's T4SS pulls non-mobilizable DNA back into the donor cell.",
        ], className="chart-subtitle"),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Transfer Routes for Non-mobilizable Plasmids",
                         className="chart-title"),
                html.P("Multiple HGT mechanisms can mobilise plasmids "
                       "that lack their own transfer genes.",
                       className="chart-subtitle"),
                dcc.Graph(figure=route_fig, config={"displayModeBar": False}),
                html.Div(className="detail-grid", children=[
                    html.Div(className="detail-item", children=[
                        html.Span("AMR at risk", className="detail-label"),
                        html.Span(f"{retro['retro_mob_amr']:,} non-mob with AMR + conj partner",
                                  className="detail-value",
                                  style={"color": COLORS["danger"]}),
                    ]),
                ]),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("Conjugative Partners in Retro-mobilization",
                         className="chart-title"),
                html.P("Which conjugative Inc groups most frequently co-exist "
                       "with non-mobilizable plasmids.",
                       className="chart-subtitle"),
                dcc.Graph(figure=inc_fig, config={"displayModeBar": False}),
            ]),
        ]),
    ])


def _build_kpc_section():
    """Build KPC transposon context section."""
    kpc = db.kpc_context_analysis()
    if not kpc or kpc["total"] == 0:
        return html.Div()

    # KPC context genes bar chart
    ctx_df = pd.DataFrame(kpc["context"])
    ctx_fig = px.bar(ctx_df, x="gene", y="cnt",
                     color_discrete_sequence=[COLORS["danger"]])
    ctx_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=380,
        xaxis_tickangle=-45, xaxis_title="", yaxis_title="Plasmids",
    )

    # KPC Inc groups bar
    inc_df = pd.DataFrame(
        sorted(kpc["inc_groups"].items(), key=lambda x: -x[1])[:12],
        columns=["Inc Group", "Plasmids"],
    )
    inc_fig = px.bar(inc_df, x="Plasmids", y="Inc Group", orientation="h",
                     color_discrete_sequence=[COLORS["accent2"]])
    inc_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=380,
        yaxis=dict(autorange="reversed"), xaxis_title="Plasmids",
    )

    # KPC mobility donut
    mob = kpc["mobility"]
    mob_fig = go.Figure(go.Pie(
        labels=[k.capitalize() for k in mob.keys()],
        values=list(mob.values()), hole=0.55,
        marker_colors=[COLORS["accent3"], COLORS["accent"], COLORS["accent4"]],
        textinfo="label+percent",
    ))
    mob_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=300, showlegend=False,
    )

    return html.Div(className="correlation-section", children=[
        html.H3(f"10. blaKPC Transposon Context ({kpc['total']:,} plasmids)",
                 className="correlation-heading"),
        html.P("NTEKPC elements — the non-Tn4401 structures carrying blaKPC — "
               "share structural features with classical transposons "
               "(inverted repeats, target site duplications, "
               "association with IS elements). Their consistent genetic context "
               "across diverse plasmid backbones and host species supports "
               "classification as true transposons rather than arbitrary "
               "genetic islands.",
               className="chart-subtitle"),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("Genes within 5 kb of blaKPC", className="chart-title"),
                html.P("blaTEM co-occurs most frequently — part of the "
                       "NTEKPC-IId element structure.",
                       className="chart-subtitle"),
                dcc.Graph(figure=ctx_fig, config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("blaKPC by Inc Group", className="chart-title"),
                dcc.Graph(figure=inc_fig, config={"displayModeBar": False}),
            ]),
        ]),
        html.Div(className="chart-grid-1", children=[
            html.Div(className="chart-card", children=[
                html.H3("blaKPC by Mobility", className="chart-title"),
                html.P(f"{mob.get('conjugative',0):,} conjugative ({mob.get('conjugative',0)/kpc['total']*100:.0f}%) — "
                       f"most KPC plasmids are self-transmissible.",
                       className="chart-subtitle"),
                dcc.Graph(figure=mob_fig, config={"displayModeBar": False}),
            ]),
        ]),
    ])


def _build_integron_ml_section():
    """Build integron ML + transposon-AMR section."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import LabelEncoder
    from sklearn.model_selection import cross_val_score
    import numpy as np

    # Build integron ML model
    rows = db.integron_ml_features()
    if not rows or len(rows) < 100:
        return html.P("Insufficient data for integron ML.", className="text-muted")

    df = pd.DataFrame(rows)
    features = []
    X_data = {}

    for col in ["mobility", "genus"]:
        le = LabelEncoder()
        X_data[col] = le.fit_transform(df[col].fillna("Unknown"))
        features.append(col)

    X_data["log_length"] = np.log1p(df["length"].fillna(0).values.astype(float))
    features.append("log_length")
    X_data["gc"] = df["gc"].fillna(0).values.astype(float)
    features.append("gc")

    for inc in ["IncFIB", "IncFII", "IncFIA", "IncI1", "IncN", "IncHI2", "IncC"]:
        X_data[inc] = df["rep_type"].fillna("").str.contains(inc).astype(int).values
        features.append(inc)

    X = pd.DataFrame(X_data)[features]
    y = df["has_integron"].values

    rf = RandomForestClassifier(n_estimators=100, max_depth=8, random_state=42)
    cv = cross_val_score(rf, X, y, cv=5, scoring="accuracy")
    rf.fit(X, y)

    imp_df = pd.DataFrame({
        "Feature": [f.replace("log_", "Log ").capitalize() for f in features],
        "Importance": rf.feature_importances_,
    }).sort_values("Importance", ascending=True)

    imp_fig = px.bar(imp_df, x="Importance", y="Feature", orientation="h",
                     color="Importance", color_continuous_scale="Viridis")
    imp_fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=400,
        coloraxis_showscale=False,
        xaxis_title="Feature Importance", yaxis_title="",
    )

    # Transposon-AMR co-occurrence heatmap
    cooccur = db.transposon_amr_cooccurrence()
    heatmap = go.Figure()
    if cooccur:
        all_amr = {}
        for fam_data in cooccur.values():
            for amr_cls, cnt in fam_data.items():
                all_amr[amr_cls] = all_amr.get(amr_cls, 0) + cnt
        top_amr = sorted(all_amr, key=all_amr.get, reverse=True)[:10]
        top_is = sorted(cooccur, key=lambda f: sum(cooccur[f].values()),
                        reverse=True)[:10]
        z = [[cooccur[is_f].get(amr_c, 0) for amr_c in top_amr] for is_f in top_is]
        heatmap = go.Figure(go.Heatmap(
            z=z, x=top_amr, y=top_is,
            colorscale="YlOrRd", texttemplate="%{z}", textfont={"size": 8},
            hovertemplate="IS: <b>%{y}</b><br>AMR: <b>%{x}</b><br>Count: %{z:,}<extra></extra>",
        ))
        heatmap.update_layout(
            template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
            margin=dict(t=10, b=10, l=10, r=10), height=400,
            xaxis=dict(tickangle=-45, tickfont=dict(size=9)),
            yaxis=dict(autorange="reversed"),
        )

    return html.Div(className="correlation-section", children=[
        html.H3("11. Integron ML & Transposon-AMR Correlations",
                 className="correlation-heading"),
        html.Div(className="chart-grid-2", children=[
            html.Div(className="chart-card", children=[
                html.H3("ML: What Predicts Integron Carriage?",
                         className="chart-title"),
                html.P(f"Random Forest (5-fold CV: {cv.mean()*100:.1f}% accuracy) "
                       f"predicting class 1 integron carriage from plasmid features.",
                       className="chart-subtitle"),
                dcc.Graph(figure=imp_fig, config={"displayModeBar": False}),
            ]),
            html.Div(className="chart-card", children=[
                html.H3("IS Families x AMR Drug Classes",
                         className="chart-title"),
                html.P("Which IS element families co-occur with which "
                       "resistance classes on the same plasmid.",
                       className="chart-subtitle"),
                dcc.Graph(figure=heatmap, config={"displayModeBar": False}),
            ]),
        ]),
    ])


@callback(
    Output("geo-inc-animation", "children"),
    Input("geo-inc-dropdown", "value"),
    prevent_initial_call=True,
)
def update_inc_animation(inc_group):
    if not inc_group:
        return html.P("Select an Inc group", className="text-muted")
    data = db.geo_temporal_inc(inc_group)
    if not data:
        return html.P(f"No geographic data for {inc_group}", className="text-muted")
    fig = _build_animated_map(data, color_scale="Tealgrn")
    return dcc.Graph(figure=fig, config={"scrollZoom": True})


@callback(
    Output("geo-mob-animation", "children"),
    Input("geo-mob-dropdown", "value"),
    prevent_initial_call=True,
)
def update_mob_animation(mob_type):
    if not mob_type:
        return html.P("Select a mobility type", className="text-muted")
    data = db.geo_temporal_mobility(mob_type)
    if not data:
        return html.P(f"No geographic data for {mob_type}", className="text-muted")
    fig = _build_animated_map(data, color_scale="Purples")
    return dcc.Graph(figure=fig, config={"scrollZoom": True})


@callback(
    Output("viewer-result", "children"),
    Input("viewer-btn", "n_clicks"),
    State("viewer-input", "value"),
    prevent_initial_call=True,
)
def visualize_plasmid(n_clicks, accession):
    if not accession or not accession.strip():
        return html.P("Please enter an accession number.", className="text-muted")

    accession = accession.strip()
    try:
        # 1. Try local SQLite database first
        plsdb_data = db.plasmid_summary(accession)
        from_ncbi = plsdb_data is None

        if plsdb_data:
            nuc = plsdb_data["nuccore"]
            tax = plsdb_data["taxonomy"]
            typ = plsdb_data["typing"]
            plasmid_feat = {
                "length": nuc.get("NUCCORE_Length", 0),
                "topology": nuc.get("NUCCORE_Topology", "circular"),
                "gc": nuc.get("NUCCORE_GC", 0),
                "description": nuc.get("NUCCORE_Description", ""),
                "organism": tax.get("TAXONOMY_species", ""),
                "accession": accession,
                "amr": [
                    {"gene": a["gene_symbol"], "product": a["gene_name"] or "",
                     "drug_class": a["drug_class"] or "", "agent": a["antimicrobial_agent"] or "",
                     "start": a["input_gene_start"] or 0, "end": a["input_gene_stop"] or 0,
                     "strand": 1 if a.get("strand_orientation", "+") == "+" else -1}
                    for a in plsdb_data["amr"]
                    if a.get("input_gene_start")
                ],
                "mob": [
                    {"element": m["element"], "biomarker": m["biomarker"] or m["element"],
                     "start": min(m["sstart"] or 0, m["send"] or 0),
                     "end": max(m["sstart"] or 0, m["send"] or 0),
                     "strand": 1 if m.get("sstrand") == "plus" else -1,
                     "identity": m["pident"]}
                    for m in plsdb_data["mob"]
                    if m.get("sstart")
                ],
                "bgc": [],
                "typing": {
                    "rep_type": typ.get("rep_type", ""),
                    "relaxase_type": typ.get("relaxase_type", ""),
                    "mpf_type": typ.get("mpf_type", ""),
                    "orit_type": typ.get("orit_type", ""),
                    "predicted_mobility": typ.get("predicted_mobility", ""),
                    "host_range": typ.get("predicted_host_range_overall_name", ""),
                    "host_range_rank": typ.get("predicted_host_range_overall_rank", ""),
                    "pmlst_scheme": typ.get("PMLST_scheme", ""),
                    "pmlst_alleles": typ.get("PMLST_alleles", ""),
                },
                "pgap_count": 0,
            }

        # 2. Always fetch GenBank for CDS + GC content (+ fallback metadata)
        gb_data = fetch_genbank_full(accession)
        if not gb_data and from_ncbi:
            return html.P(
                f"Accession '{accession}' not found in PLSDB or NCBI GenBank.",
                className="text-danger",
            )

        gb_features = gb_data["features"] if gb_data else []
        gc_profile = gb_data.get("gc_profile", []) if gb_data else []

        if from_ncbi:
            gb_meta = gb_data["metadata"]
            plasmid_feat = {
                "length": gb_meta.get("length", 0),
                "topology": gb_meta.get("topology", "circular"),
                "gc": gb_meta.get("gc_overall", 0),
                "description": gb_meta.get("description", ""),
                "organism": gb_meta.get("organism", ""),
                "accession": gb_meta.get("accession", accession),
                "amr": [], "mob": [], "bgc": [],
                "typing": {}, "pgap_count": 0,
            }

        # Use GC from GenBank if DB didn't provide it
        if plasmid_feat.get("gc", 0) == 0 and gb_data:
            plasmid_feat["gc"] = gb_data["metadata"].get("gc_overall", 0)

        # Classify all GenBank CDS genes into functional categories
        gene_cats = {}  # category -> list of gene names
        if gb_features:
            for f in gb_features:
                gene = f.get("gene") or f.get("locus_tag") or ""
                cat = classify_gene(gene, f.get("product", ""))
                if cat != "CDS":
                    gene_cats.setdefault(cat, []).append(gene)

        # Build the interactive circular map (Plotly figure)
        fig = make_plasmid_map(plasmid_feat, gb_features, gc_profile)

        # Build typing info card
        typing = plasmid_feat.get("typing", {})
        typing_items = []
        typing_fields = [
            ("Inc Group / Replicon Type", typing.get("rep_type")),
            ("Relaxase Type", typing.get("relaxase_type")),
            ("MPF Type", typing.get("mpf_type")),
            ("OriT Type", typing.get("orit_type")),
            ("Predicted Mobility", typing.get("predicted_mobility")),
            ("Host Range", f"{typing.get('host_range', '')} ({typing.get('host_range_rank', '')})"),
            ("pMLST", typing.get("pmlst_alleles")),
        ]
        for label, value in typing_fields:
            if value and str(value).strip() and str(value).strip() != "()":
                typing_items.append(
                    html.Div(className="detail-item", children=[
                        html.Span(label, className="detail-label"),
                        html.Span(str(value), className="detail-value"),
                    ])
                )

        # Feature summary badges
        summary_items = []
        if gb_features:
            summary_items.append(f"{len(gb_features)} CDS genes")

        # AMR: combine PLSDB + GenBank detected
        plsdb_amr = plasmid_feat.get("amr", [])
        gb_amr = gene_cats.get("AMR", [])
        if plsdb_amr:
            summary_items.append(f"{len(plsdb_amr)} AMR (PLSDB)")
        if gb_amr:
            summary_items.append(f"{len(gb_amr)} AMR: {', '.join(gb_amr)}")

        # MOB: combine PLSDB + GenBank detected
        plsdb_mob = plasmid_feat.get("mob", [])
        gb_mob = gene_cats.get("MOB", [])
        if plsdb_mob:
            summary_items.append(f"{len(plsdb_mob)} MOB (PLSDB)")
        if gb_mob:
            summary_items.append(f"{len(gb_mob)} MOB: {', '.join(gb_mob)}")

        # Toxin-Antitoxin
        gb_ta = gene_cats.get("TA", [])
        if gb_ta:
            summary_items.append(f"{len(gb_ta)} TA: {', '.join(gb_ta)}")

        # Virulence
        gb_vir = gene_cats.get("VIR", [])
        if gb_vir:
            summary_items.append(f"{len(gb_vir)} Virulence: {', '.join(gb_vir)}")

        # QAC
        gb_qac = gene_cats.get("QAC", [])
        if gb_qac:
            summary_items.append(f"{len(gb_qac)} QAC: {', '.join(gb_qac)}")

        # Heavy metal
        gb_metal = gene_cats.get("Metal", [])
        if gb_metal:
            summary_items.append(f"{len(gb_metal)} Metal: {', '.join(gb_metal)}")

        # BGC
        plsdb_bgc = plasmid_feat.get("bgc", [])
        if plsdb_bgc:
            summary_items.append(f"{len(plsdb_bgc)} BGC regions")

        if plasmid_feat.get("pgap_count"):
            summary_items.append(f"{plasmid_feat['pgap_count']} annotations (PGAP)")

        # NCBI fallback notice
        ncbi_notice = []
        if from_ncbi:
            ncbi_notice = [html.P(
                "Not found in PLSDB. Showing GenBank data from NCBI.",
                className="text-muted",
                style={"fontStyle": "italic", "marginBottom": "8px"},
            )]

        # Link
        if from_ncbi:
            ext_link = html.A(
                f"View {accession} in NCBI",
                href=f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                target="_blank", className="btn-secondary",
                style={"marginTop": "15px", "display": "inline-block"},
            )
        else:
            ext_link = html.A(
                f"View {accession} in PLSDB",
                href=f"https://ccb-microbe.cs.uni-saarland.de/plsdb2025/plasmid/{accession}",
                target="_blank", className="btn-secondary",
                style={"marginTop": "15px", "display": "inline-block"},
            )

        return html.Div(ncbi_notice + [
            html.Div(className="feature-badges", children=[
                html.Span(s, className="badge") for s in summary_items
            ]),
            # Interactive plasmid map (Plotly with zoom + hover)
            dcc.Graph(
                figure=fig,
                config={"scrollZoom": True, "displayModeBar": True},
                style={"maxWidth": "850px", "margin": "0 auto"},
            ),
            # Typing info
            html.H4("MOB Typing & Incompatibility", className="chart-title",
                     style={"marginTop": "20px"}),
            html.Div(className="detail-grid", children=typing_items) if typing_items
            else html.P("No MOB typing data available.", className="text-muted"),
            ext_link,
        ])

    except Exception as e:
        logger.error(f"Viewer error: {e}", exc_info=True)
        return html.P(f"Error: {str(e)}", className="text-danger")


@callback(
    Output("lookup-result", "children"),
    Input("lookup-btn", "n_clicks"),
    State("lookup-input", "value"),
    prevent_initial_call=True,
)
def lookup_plasmid(n_clicks, accession):
    if not accession or not accession.strip():
        return html.P("Please enter an accession number.", className="text-muted")

    accession = accession.strip()
    try:
        data = db.plasmid_summary(accession)
        if not data:
            return html.P(f"Accession '{accession}' not found in PLSDB.",
                          className="text-danger")

        nuc = data["nuccore"]
        tax = data["taxonomy"]

        detail_fields = [
            ("Accession", nuc.get("NUCCORE_ACC")),
            ("Description", nuc.get("NUCCORE_Description")),
            ("Topology", nuc.get("NUCCORE_Topology")),
            ("Length (bp)", f"{nuc.get('NUCCORE_Length', 0):,}"),
            ("GC Content", f"{(nuc.get('NUCCORE_GC', 0) or 0) * 100:.1f}%"),
            ("Source", nuc.get("NUCCORE_Source")),
            ("Create Date", nuc.get("NUCCORE_CreateDate")),
            ("Species", tax.get("TAXONOMY_species")),
            ("Genus", tax.get("TAXONOMY_genus")),
            ("Family", tax.get("TAXONOMY_family")),
            ("Phylum", tax.get("TAXONOMY_phylum")),
        ]

        detail_items = []
        for label, value in detail_fields:
            if value and str(value).strip():
                detail_items.append(
                    html.Div(className="detail-item", children=[
                        html.Span(label, className="detail-label"),
                        html.Span(str(value), className="detail-value"),
                    ])
                )

        if detail_items:
            cards.append(html.Div(className="detail-grid", children=detail_items))

        # Link to PLSDB page
        cards.append(
            html.A(
                f"View {accession} in PLSDB",
                href=f"https://ccb-microbe.cs.uni-saarland.de/plsdb2025/plasmid/{accession}",
                target="_blank", className="btn-secondary",
                style={"marginTop": "15px", "display": "inline-block"},
            )
        )

        return html.Div(cards)

    except Exception as e:
        return html.P(f"Error: {str(e)}", className="text-danger")


# ---------------------------------------------------------------------------
# Pre-compute ALL expensive charts (runs once at startup, tabs load instantly)
# ---------------------------------------------------------------------------
logger.info("Pre-computing all charts...")
_shap_fig, _shap_acc = make_xgboost_shap_chart()
logger.info(f"  SHAP done ({_shap_acc}%)")

# Cache all chart figures as module-level variables
_CACHE = {}
_cache_fns = [
    ("rf", make_feature_importance_chart),
    ("matched", make_matched_comparison),
    ("rarefaction", make_rarefaction_chart),
    ("amr_trends", make_temporal_trends_chart),
    ("inc_trends", make_inc_trends_chart),
    ("amr_cooccur", make_amr_cooccurrence_chart),
    ("umap_mob", lambda: make_umap_chart("mobility")),
    ("umap_genus", lambda: make_umap_chart("genus")),
    ("simpson", make_simpson_paradox_table),
    ("global_map", make_global_map),
    ("country_mob", make_country_mobility_chart),
    ("country_inc", make_country_inc_chart),
    ("temporal_anim", make_temporal_animation),
    ("host_dist", make_host_distribution_chart),
    ("host_mob", make_host_mobility_chart),
    ("host_inc", make_host_inc_heatmap),
    ("host_amr", make_host_amr_heatmap),
    ("is_family", make_is_family_chart),
    ("is_mob", make_is_mobility_chart),
    ("is_inc", make_is_inc_heatmap),
    ("integron_sec", _build_integron_section),
    ("comob_sec", _build_comob_section),
    ("retromob_sec", _build_retromob_section),
    ("kpc_sec", _build_kpc_section),
    ("integron_ml_sec", _build_integron_ml_section),
]
for _name, _fn in _cache_fns:
    logger.info(f"  {_name}...")
    _CACHE[_name] = _fn()

# Override all heavy functions to return cached results
make_feature_importance_chart = lambda: _CACHE["rf"]
make_matched_comparison = lambda: _CACHE["matched"]
make_rarefaction_chart = lambda: _CACHE["rarefaction"]
make_temporal_trends_chart = lambda: _CACHE["amr_trends"]
make_inc_trends_chart = lambda: _CACHE["inc_trends"]
make_amr_cooccurrence_chart = lambda: _CACHE["amr_cooccur"]
make_umap_chart = lambda c="mobility": _CACHE["umap_mob"] if c == "mobility" else _CACHE["umap_genus"]
make_simpson_paradox_table = lambda: _CACHE["simpson"]
make_global_map = lambda: _CACHE["global_map"]
make_country_mobility_chart = lambda: _CACHE["country_mob"]
make_country_inc_chart = lambda: _CACHE["country_inc"]
make_temporal_animation = lambda: _CACHE["temporal_anim"]
make_host_distribution_chart = lambda: _CACHE["host_dist"]
make_host_mobility_chart = lambda: _CACHE["host_mob"]
make_host_inc_heatmap = lambda: _CACHE["host_inc"]
make_host_amr_heatmap = lambda: _CACHE["host_amr"]
make_is_family_chart = lambda: _CACHE["is_family"]
make_is_mobility_chart = lambda: _CACHE["is_mob"]
make_is_inc_heatmap = lambda: _CACHE["is_inc"]
_build_integron_section = lambda: _CACHE["integron_sec"]
_build_comob_section = lambda: _CACHE["comob_sec"]
_build_retromob_section = lambda: _CACHE["retromob_sec"]
_build_kpc_section = lambda: _CACHE["kpc_sec"]
_build_integron_ml_section = lambda: _CACHE["integron_ml_sec"]

logger.info("All charts pre-computed! Tabs will load instantly.")

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8050))
    debug = os.environ.get("DASH_DEBUG", "true").lower() == "true"
    app.run(debug=debug, host="0.0.0.0", port=port)
