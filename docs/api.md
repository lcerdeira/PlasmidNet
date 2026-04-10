# REST API

PlasmidNet provides a REST API for programmatic access to the plasmid database. All endpoints return JSON.

**Base URL**: `https://plasmidnet-7302570f4818.herokuapp.com/api/`

## Endpoints

### GET /api/

API documentation and endpoint list.

```bash
curl https://plasmidnet-7302570f4818.herokuapp.com/api/
```

```json
{
  "name": "PlasmidNet API",
  "version": "1.0",
  "total_plasmids": 72556,
  "endpoints": {
    "/api/plasmid/<accession>": "Get plasmid details",
    "/api/search?q=<query>": "Search plasmids",
    "/api/amr?gene=<gene>": "AMR gene search",
    "/api/amr/classes": "Drug class counts",
    "/api/typing?inc=<inc>": "Inc group search",
    "/api/stats": "Database statistics"
  }
}
```

### GET /api/plasmid/{accession}

Get full details for a single plasmid including nuccore, typing, and AMR data.

```bash
curl https://plasmidnet-7302570f4818.herokuapp.com/api/plasmid/NZ_CP031107.1
```

**Response**: JSON with `nuccore`, `taxonomy`, `typing`, `amr`, `mob` fields.

### GET /api/search

Search plasmids by description, organism, or accession.

| Parameter | Type | Description |
|-----------|------|-------------|
| `q` | string | Search query (min 2 chars) |
| `limit` | int | Max results (default 100, max 500) |

```bash
# Search for KPC plasmids
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/search?q=KPC&limit=10"

# Search by organism
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/search?q=Klebsiella"
```

### GET /api/amr

Get plasmids carrying a specific AMR gene or drug class.

| Parameter | Type | Description |
|-----------|------|-------------|
| `gene` | string | AMR gene symbol (e.g., `blaTEM`, `sul1`) |
| `class` | string | Drug class (e.g., `BETA-LACTAM`, `MERCURY`) |
| `limit` | int | Max results (default 100) |

```bash
# Find plasmids with blaKPC
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/amr?gene=blaKPC"

# Find mercury resistance
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/amr?class=MERCURY"
```

### GET /api/amr/classes

List all AMR drug classes with plasmid counts.

```bash
curl https://plasmidnet-7302570f4818.herokuapp.com/api/amr/classes
```

### GET /api/amr/genes

Top AMR genes with counts.

```bash
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/amr/genes?limit=10"
```

### GET /api/typing

Get plasmids by Inc group or mobility type.

| Parameter | Type | Description |
|-----------|------|-------------|
| `inc` | string | Inc group (e.g., `IncFIB`, `IncFII`) |
| `mobility` | string | `conjugative`, `mobilizable`, or `non-mobilizable` |
| `limit` | int | Max results (default 100) |

```bash
# Find IncFIB plasmids
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/typing?inc=IncFIB&limit=5"

# Find conjugative plasmids
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/typing?mobility=conjugative&limit=5"
```

### GET /api/typing/mobility

Get mobility distribution.

```bash
curl https://plasmidnet-7302570f4818.herokuapp.com/api/typing/mobility
```

```json
{
  "conjugative": 22588,
  "mobilizable": 21152,
  "non-mobilizable": 28816
}
```

### GET /api/typing/inc

Top Inc groups with counts.

```bash
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/typing/inc?limit=10"
```

### GET /api/taxonomy

Get plasmids by genus or species.

| Parameter | Type | Description |
|-----------|------|-------------|
| `genus` | string | Genus name (e.g., `Escherichia`) |
| `species` | string | Species name (e.g., `Escherichia coli`) |
| `limit` | int | Max results (default 100) |

```bash
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/taxonomy?genus=Klebsiella&limit=5"
```

### GET /api/taxonomy/genera

Top genera with counts.

```bash
curl "https://plasmidnet-7302570f4818.herokuapp.com/api/taxonomy/genera?limit=10"
```

### GET /api/stats

Database summary statistics.

```bash
curl https://plasmidnet-7302570f4818.herokuapp.com/api/stats
```

## Python Example

```python
import requests

BASE = "https://plasmidnet-7302570f4818.herokuapp.com/api"

# Get all KPC plasmids
kpc = requests.get(f"{BASE}/amr", params={"gene": "blaKPC"}).json()
print(f"Found {kpc['count']} plasmids with blaKPC")

# Get details for a specific plasmid
p = requests.get(f"{BASE}/plasmid/NZ_CP031107.1").json()
print(f"Length: {p['nuccore']['NUCCORE_Length']} bp")
print(f"Mobility: {p['typing']['predicted_mobility']}")
print(f"AMR genes: {len(p['amr'])}")

# Search
results = requests.get(f"{BASE}/search", params={"q": "carbapenem"}).json()
for r in results["results"][:5]:
    print(f"  {r['NUCCORE_ACC']}: {r['NUCCORE_Description'][:60]}")
```

## R Example

```r
library(httr)
library(jsonlite)

base <- "https://plasmidnet-7302570f4818.herokuapp.com/api"

# Get mobility distribution
mob <- fromJSON(content(GET(paste0(base, "/typing/mobility")), "text"))
print(mob)

# Search for Salmonella plasmids
sal <- fromJSON(content(GET(paste0(base, "/taxonomy?genus=Salmonella&limit=5")), "text"))
print(sal$results)
```

## Rate Limits

No authentication required. Please limit requests to 10/second to avoid overloading the Heroku dyno.
