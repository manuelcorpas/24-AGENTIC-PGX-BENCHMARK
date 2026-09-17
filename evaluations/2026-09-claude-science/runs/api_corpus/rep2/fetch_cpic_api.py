#!/usr/bin/env python3
"""Snapshot the CPIC data API (api.cpicpgx.org/v1) tables used by pgx_pipeline.py.

Writes cpic_api_snapshot/<table>.json, cpic_api_snapshot/diplotype_by_gene.json and
cpic_api_snapshot/_fetch_meta.json (fetch timestamp + sha256 of every file).
"""
import json, os, hashlib, datetime, requests

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = os.path.join(HERE, "cpic_api_snapshot")
BASE = "https://api.cpicpgx.org/v1/"
TABLES = ["gene", "pair", "drug", "guideline", "allele", "allele_definition",
          "allele_location_value", "sequence_location", "gene_result",
          "gene_result_lookup", "recommendation", "term", "change_log"]

def fetch(table, **params):
    p = {"limit": "50000"}
    p.update({k: str(v) for k, v in params.items()})
    r = requests.get(BASE + table, params=p, timeout=180)
    r.raise_for_status()
    return r.json()

def fetch_all(table, **params):
    rows, off = [], 0
    while True:
        chunk = fetch(table, offset=off, **params)
        rows += chunk
        off += len(chunk)
        if len(chunk) < 50000:
            return rows

def main():
    os.makedirs(SNAP, exist_ok=True)
    genes = sorted({c["gene"] for c in json.load(open(os.path.join(HERE, "cases_blinded.json")))})
    for t in TABLES:
        json.dump(fetch_all(t), open(os.path.join(SNAP, t + ".json"), "w"))
    # diplotype is large: fetch per gene with pagination
    dip = {g: fetch_all("diplotype", genesymbol="eq." + g) for g in genes}
    json.dump(dip, open(os.path.join(SNAP, "diplotype_by_gene.json"), "w"))

    files = sorted(f for f in os.listdir(SNAP) if f.endswith(".json") and not f.startswith("_"))
    meta = {
        "source": "CPIC data API",
        "base_url": BASE,
        "fetched_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "study_genes": genes,
        "files": {f: {"sha256": hashlib.sha256(open(os.path.join(SNAP, f), "rb").read()).hexdigest(),
                      "bytes": os.path.getsize(os.path.join(SNAP, f))} for f in files},
    }
    cl = json.load(open(os.path.join(SNAP, "change_log.json")))
    dates = sorted({c.get("date") for c in cl if c.get("date")})
    meta["change_log_latest_date"] = dates[-1] if dates else None
    meta["table_row_counts"] = {f[:-5]: len(json.load(open(os.path.join(SNAP, f)))) for f in files}
    json.dump(meta, open(os.path.join(SNAP, "_fetch_meta.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != "files"}, indent=1))

if __name__ == "__main__":
    main()
