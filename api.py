"""
PlasmidNet REST API — programmatic access to the SQLite database.
Mounted at /api/ in the Dash app.
"""

from flask import Blueprint, jsonify, request
import db

api_bp = Blueprint("api", __name__, url_prefix="/api")


@api_bp.route("/")
def api_docs():
    """API documentation."""
    return jsonify({
        "name": "PlasmidNet API",
        "version": "1.0",
        "total_plasmids": db.scalar("SELECT COUNT(*) FROM nuccore"),
        "endpoints": {
            "/api/plasmid/<accession>": "Get plasmid details (nuccore + typing + AMR)",
            "/api/search?q=<query>": "Search plasmids by description, organism, or accession",
            "/api/amr?gene=<gene>": "Get plasmids carrying a specific AMR gene",
            "/api/amr/classes": "List all AMR drug classes with counts",
            "/api/typing?inc=<inc_group>": "Get plasmids by Inc group",
            "/api/typing/mobility": "Mobility distribution",
            "/api/taxonomy?genus=<genus>": "Get plasmids by genus",
            "/api/stats": "Database summary statistics",
        },
    })


@api_bp.route("/plasmid/<accession>")
def get_plasmid(accession):
    """Get full plasmid details."""
    data = db.plasmid_summary(accession)
    if not data:
        return jsonify({"error": f"Accession {accession} not found"}), 404
    return jsonify(data)


@api_bp.route("/search")
def search_plasmids():
    """Search plasmids by description, organism, or accession."""
    q = request.args.get("q", "").strip()
    limit = min(int(request.args.get("limit", 100)), 500)
    if not q or len(q) < 2:
        return jsonify({"error": "Query must be at least 2 characters"}), 400

    rows = db.q("""
        SELECT n.NUCCORE_ACC, n.NUCCORE_Description, n.NUCCORE_Length,
               n.NUCCORE_Topology, n.NUCCORE_GC,
               t.TAXONOMY_genus, t.TAXONOMY_species
        FROM nuccore n
        LEFT JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
        WHERE n.NUCCORE_ACC LIKE ? OR n.NUCCORE_Description LIKE ?
           OR t.TAXONOMY_genus LIKE ? OR t.TAXONOMY_species LIKE ?
        LIMIT ?
    """, (f"%{q}%", f"%{q}%", f"%{q}%", f"%{q}%", limit))

    return jsonify({"query": q, "count": len(rows), "results": rows})


@api_bp.route("/amr")
def get_amr():
    """Get plasmids carrying a specific AMR gene or drug class."""
    gene = request.args.get("gene", "").strip()
    drug_class = request.args.get("class", "").strip()
    limit = min(int(request.args.get("limit", 100)), 500)

    if gene:
        rows = db.q("""
            SELECT NUCCORE_ACC, gene_symbol, gene_name, drug_class,
                   input_gene_start, input_gene_stop
            FROM amr WHERE gene_symbol LIKE ?
            LIMIT ?
        """, (f"%{gene}%", limit))
    elif drug_class:
        rows = db.q("""
            SELECT NUCCORE_ACC, gene_symbol, gene_name, drug_class
            FROM amr WHERE UPPER(drug_class) LIKE ?
            LIMIT ?
        """, (f"%{drug_class.upper()}%", limit))
    else:
        return jsonify({"error": "Provide ?gene= or ?class= parameter"}), 400

    return jsonify({"count": len(rows), "results": rows})


@api_bp.route("/amr/classes")
def amr_classes():
    """List all AMR drug classes with counts."""
    return jsonify(db.amr_drug_class_counts())


@api_bp.route("/amr/genes")
def amr_genes():
    """Top AMR genes with counts."""
    limit = min(int(request.args.get("limit", 30)), 100)
    return jsonify(db.amr_gene_counts(limit))


@api_bp.route("/typing")
def get_typing():
    """Get plasmids by Inc group."""
    inc = request.args.get("inc", "").strip()
    mobility = request.args.get("mobility", "").strip()
    limit = min(int(request.args.get("limit", 100)), 500)

    if inc:
        rows = db.q("""
            SELECT NUCCORE_ACC, rep_type, predicted_mobility,
                   relaxase_type, mpf_type
            FROM typing WHERE rep_type LIKE ?
            LIMIT ?
        """, (f"%{inc}%", limit))
    elif mobility:
        rows = db.q("""
            SELECT NUCCORE_ACC, rep_type, predicted_mobility,
                   relaxase_type, mpf_type
            FROM typing WHERE predicted_mobility = ?
            LIMIT ?
        """, (mobility, limit))
    else:
        return jsonify({"error": "Provide ?inc= or ?mobility= parameter"}), 400

    return jsonify({"count": len(rows), "results": rows})


@api_bp.route("/typing/mobility")
def mobility_dist():
    """Mobility distribution."""
    return jsonify(db.mobility_distribution())


@api_bp.route("/typing/inc")
def inc_groups():
    """Inc group counts."""
    limit = min(int(request.args.get("limit", 20)), 50)
    return jsonify(db.inc_group_counts(limit))


@api_bp.route("/taxonomy")
def get_taxonomy():
    """Get plasmids by genus or species."""
    genus = request.args.get("genus", "").strip()
    species = request.args.get("species", "").strip()
    limit = min(int(request.args.get("limit", 100)), 500)

    where = []
    params = []
    if genus:
        where.append("t.TAXONOMY_genus LIKE ?")
        params.append(f"%{genus}%")
    if species:
        where.append("t.TAXONOMY_species LIKE ?")
        params.append(f"%{species}%")

    if not where:
        return jsonify({"error": "Provide ?genus= or ?species= parameter"}), 400

    rows = db.q(f"""
        SELECT n.NUCCORE_ACC, n.NUCCORE_Description, n.NUCCORE_Length,
               t.TAXONOMY_genus, t.TAXONOMY_species
        FROM nuccore n
        JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
        WHERE {' AND '.join(where)}
        LIMIT ?
    """, tuple(params) + (limit,))

    return jsonify({"count": len(rows), "results": rows})


@api_bp.route("/taxonomy/genera")
def top_genera():
    """Top genera with counts."""
    limit = min(int(request.args.get("limit", 20)), 50)
    return jsonify(db.top_genera(limit))


@api_bp.route("/stats")
def stats():
    """Database summary statistics."""
    overview = db.overview_stats()
    overview["mobility"] = db.mobility_distribution()
    return jsonify(overview)
