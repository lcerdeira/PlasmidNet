"""
PlasmidNet Dashboard - Interactive explorer for the PLSDB 2025 plasmid database.
Built with Dash + Plotly. Deployable on Heroku / AWS free tier.
"""

import os
import logging

import dash
from dash import dcc, html, dash_table, callback, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import requests

from data_loader import load_all_data, fetch_summary, API_BASE

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
DATA = load_all_data(use_fallback=True)
overview = DATA["overview"]
taxonomy = DATA["taxonomy"]
amr = DATA["amr"]
temporal = DATA.get("temporal", {})
gc_dist = DATA.get("gc_distribution", {})
length_dist = DATA.get("length_distribution", {})

# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------
app = dash.Dash(
    __name__,
    title="PlasmidNet - PLSDB 2025 Dashboard",
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    suppress_callback_exceptions=True,
)
server = app.server  # for gunicorn

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
COLORS = {
    "bg": "#0f172a",
    "card": "#1e293b",
    "card_border": "#334155",
    "accent": "#38bdf8",
    "accent2": "#818cf8",
    "accent3": "#34d399",
    "accent4": "#fb923c",
    "accent5": "#f472b6",
    "text": "#f1f5f9",
    "text_muted": "#94a3b8",
    "danger": "#ef4444",
}

PLOTLY_TEMPLATE = "plotly_dark"
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
# Layout
# ---------------------------------------------------------------------------

HEADER = html.Div(className="header", children=[
    html.Div(className="header-content", children=[
        html.Div([
            html.H1("PlasmidNet", className="logo"),
            html.Span("PLSDB 2025 Dashboard", className="subtitle"),
        ]),
        html.Div(className="header-links", children=[
            html.A("PLSDB Database", href="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/",
                    target="_blank", className="header-link"),
            html.A("API Docs", href="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api/",
                    target="_blank", className="header-link"),
            html.A("GitHub", href="https://github.com/CCB-SB/plsdbapi",
                    target="_blank", className="header-link"),
        ]),
    ]),
])

STATS_ROW = html.Div(className="stats-row", children=[
    stat_card("Total Plasmids", overview.get("total", 72360), "🧬", COLORS["accent"]),
    stat_card("Annotations", overview.get("total_annotations", 6027698), "📋", COLORS["accent2"]),
    stat_card("Virulence Factors", overview.get("total_virulence_factors", 250691), "⚠", COLORS["accent4"]),
    stat_card("Biosynthetic Gene Clusters", overview.get("total_bgcs", 12796), "🔬", COLORS["accent3"]),
])

TABS = html.Div(className="tabs-container", children=[
    dcc.Tabs(
        id="main-tabs",
        value="overview",
        children=[
            dcc.Tab(label="Overview", value="overview", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="Taxonomy", value="taxonomy", className="tab", selected_className="tab--selected"),
            dcc.Tab(label="AMR Analysis", value="amr", className="tab", selected_className="tab--selected"),
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
        "Data from ",
        html.A("PLSDB 2025", href="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/",
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
