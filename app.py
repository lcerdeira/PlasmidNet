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
import requests

from data_loader import (
    load_all_data, fetch_summary, fetch_genbank_features,
    fetch_genbank_full, extract_plasmid_features, API_BASE,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
DATA = load_all_data(use_fallback=True)
overview = DATA["overview"]
taxonomy = DATA["taxonomy"]
amr = DATA["amr"]
mob_typing = DATA.get("mob_typing", {})
correlations = DATA.get("correlations", {})
temporal = DATA.get("temporal", {})
gc_dist = DATA.get("gc_distribution", {})
length_dist = DATA.get("length_distribution", {})

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
    if not amr:
        return go.Figure()
    df = pd.DataFrame(
        sorted(amr.items(), key=lambda x: x[1], reverse=True),
        columns=["Gene", "Plasmid Count"],
    )
    # Categorize by drug class
    drug_class_map = {
        "blaTEM": "Beta-lactams", "blaSHV": "Beta-lactams", "blaCTX-M": "Beta-lactams",
        "blaOXA": "Beta-lactams", "blaKPC": "Carbapenems", "blaNDM": "Carbapenems",
        "mcr-1": "Colistin", "vanA": "Glycopeptides", "vanB": "Glycopeptides",
        "mecA": "Methicillin", "ermB": "Macrolides", "ermC": "Macrolides",
        "tetA": "Tetracyclines", "tetB": "Tetracyclines", "tetM": "Tetracyclines",
        "sul1": "Sulfonamides", "sul2": "Sulfonamides",
        "aph(3')-Ia": "Aminoglycosides", "aac(6')-Ib": "Aminoglycosides",
        "qnrS": "Quinolones", "qnrB": "Quinolones",
        "dfrA1": "Trimethoprim", "dfrA12": "Trimethoprim",
        "catA1": "Chloramphenicol", "floR": "Chloramphenicol",
        "fosA": "Fosfomycin", "aadA1": "Aminoglycosides",
        "strA": "Aminoglycosides", "strB": "Aminoglycosides",
    }
    df["Drug Class"] = df["Gene"].map(drug_class_map).fillna("Other")
    fig = px.bar(
        df, x="Gene", y="Plasmid Count", color="Drug Class",
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=10, l=10, r=10), height=400,
        xaxis_title="", yaxis_title="Number of Plasmids",
        xaxis_tickangle=-45, legend=dict(
            orientation="h", yanchor="bottom", y=1.02,
            xanchor="right", x=1, font_size=10,
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
    """Aggregate AMR genes by drug class."""
    if not amr:
        return go.Figure()
    drug_class_map = {
        "blaTEM": "Beta-lactams", "blaSHV": "Beta-lactams", "blaCTX-M": "Beta-lactams",
        "blaOXA": "Beta-lactams", "blaKPC": "Carbapenems", "blaNDM": "Carbapenems",
        "mcr-1": "Colistin", "vanA": "Glycopeptides", "vanB": "Glycopeptides",
        "mecA": "Methicillin", "ermB": "Macrolides", "ermC": "Macrolides",
        "tetA": "Tetracyclines", "tetB": "Tetracyclines", "tetM": "Tetracyclines",
        "sul1": "Sulfonamides", "sul2": "Sulfonamides",
        "aph(3')-Ia": "Aminoglycosides", "aac(6')-Ib": "Aminoglycosides",
        "qnrS": "Quinolones", "qnrB": "Quinolones",
        "dfrA1": "Trimethoprim", "dfrA12": "Trimethoprim",
        "catA1": "Chloramphenicol", "floR": "Chloramphenicol",
        "fosA": "Fosfomycin", "aadA1": "Aminoglycosides",
        "strA": "Aminoglycosides", "strB": "Aminoglycosides",
    }
    class_totals = {}
    for gene, count in amr.items():
        cls = drug_class_map.get(gene, "Other")
        class_totals[cls] = class_totals.get(cls, 0) + count

    df = pd.DataFrame(
        sorted(class_totals.items(), key=lambda x: x[1], reverse=True),
        columns=["Drug Class", "Plasmid Count"],
    )
    fig = px.bar(
        df, x="Plasmid Count", y="Drug Class", orientation="h",
        color="Drug Class", color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        template=PLOTLY_TEMPLATE, paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font_color=COLORS["text"],
        margin=dict(t=10, b=30, l=10, r=10), height=400,
        yaxis=dict(autorange="reversed"), showlegend=False,
        xaxis_title="Total Plasmids with Resistance", yaxis_title="",
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
# Layout
# ---------------------------------------------------------------------------

HEADER = html.Div(className="header", children=[
    html.Div(className="header-content", children=[
        html.Div([
            html.Img(src="/assets/logo.png", className="logo-img"),
        ]),
        html.Div(className="header-links", children=[
            html.A("Database", href="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/",
                    target="_blank", className="header-link"),
            html.A("Documentation", href="https://github.com/lcerdeira/plasmidnet/wiki",
                    target="_blank", className="header-link"),
            html.A("GitHub", href="https://github.com/lcerdeira/plasmidnet",
                    target="_blank", className="header-link"),
        ]),
    ]),
])

STATS_ROW = html.Div(className="stats-row", children=[
    stat_card("Total Plasmids", overview.get("total", 72360), COLORS["accent"]),
    stat_card("Annotations", overview.get("total_annotations", 6027698), COLORS["accent2"]),
    stat_card("Virulence Factors", overview.get("total_virulence_factors", 250691), COLORS["accent4"]),
    stat_card("Biosynthetic Gene Clusters", overview.get("total_bgcs", 12796), COLORS["accent3"]),
])

TABS = html.Div(className="tabs-container", children=[
    dcc.Tabs(
        id="main-tabs",
        value="overview",
        children=[
            dcc.Tab(label="Overview", value="overview", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Taxonomy", value="taxonomy", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="AMR Analysis", value="amr", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Plasmid Viewer", value="viewer", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Inc Groups & Mobility", value="mobility", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Correlations", value="correlations", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Plasmid Lookup", value="lookup", className="tab", selected_className="tab--selected"),
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

        # --- 3. Phage & Mobile Elements ---
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
                    html.P("Conjugative plasmids carry more phage and integrase "
                           "elements on average",
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
                    html.P("Toxin-antitoxin system components",
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

app.layout = html.Div(className="app-container", children=[
    HEADER,
    html.Div(className="main-content", children=[
        STATS_ROW,
        TABS,
        html.Div(id="tab-content"),
    ]),
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
    elif tab == "lookup":
        return lookup_tab()
    return overview_tab()


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

    params = {}
    if genus:
        params["TAXONOMY_genus"] = genus
    if species:
        params["TAXONOMY_species"] = species

    try:
        resp = requests.get(f"{API_BASE}/filter_taxonomy", params=params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, list):
                count = len(data)
            elif isinstance(data, dict):
                # Try to extract count
                for key in ["results", "data", "accessions"]:
                    if key in data:
                        count = len(data[key])
                        break
                else:
                    count = len(data)
            else:
                count = 0

            search_term = species or genus
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
            return html.P(f"No results found (status {resp.status_code})", className="text-muted")
    except Exception as e:
        return html.P(f"Error querying API: {str(e)}", className="text-danger")


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
        # 1. Try PLSDB first
        data = fetch_summary(accession)
        from_ncbi = False

        # Fix: PLSDB returns {"count":0} for missing accessions, not "notfound"
        plsdb_found = (data and "Metadata_annotations" in data)

        if plsdb_found:
            plasmid_feat = extract_plasmid_features(data)
        else:
            from_ncbi = True

        # 2. Always fetch GenBank for CDS + GC content + fallback metadata
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

        # Use GC from GenBank if PLSDB didn't provide it
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
        data = fetch_summary(accession)
        if not data or data.get("label") == "notfound":
            return html.P(f"Accession '{accession}' not found in PLSDB.",
                          className="text-danger")

        # Build detail cards from the response
        cards = []

        # Extract metadata - the response structure varies
        meta = data.get("Metadata_annotations", data)

        detail_fields = [
            ("Accession", "NUCCORE_ACC"),
            ("Description", "NUCCORE_Description"),
            ("Topology", "NUCCORE_Topology"),
            ("Length (bp)", "NUCCORE_Length"),
            ("GC Content (%)", "NUCCORE_GC"),
            ("Source", "NUCCORE_Source"),
            ("Create Date", "NUCCORE_CreateDate"),
            ("Organism", "TAXONOMY_species"),
            ("Genus", "TAXONOMY_genus"),
            ("Family", "TAXONOMY_family"),
            ("Phylum", "TAXONOMY_phylum"),
        ]

        detail_items = []
        for label, key in detail_fields:
            value = _deep_get(meta, key)
            if value is not None and str(value).strip():
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


def _deep_get(d, key):
    """Recursively search for a key in nested dicts."""
    if not isinstance(d, dict):
        return None
    if key in d:
        return d[key]
    for v in d.values():
        if isinstance(v, dict):
            result = _deep_get(v, key)
            if result is not None:
                return result
    return None


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8050))
    debug = os.environ.get("DASH_DEBUG", "true").lower() == "true"
    app.run(debug=debug, host="0.0.0.0", port=port)
