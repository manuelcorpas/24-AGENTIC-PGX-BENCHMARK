#!/usr/bin/env python3
"""
Download the CPIC data API tables this study depends on into ./cpic_api_raw/.

Run once; pgx_pipeline.py reads the cached JSON so the analysis is reproducible
without network access.  api.cpicpgx.org is a PostgREST endpoint: tables are
paths and filters are query parameters (e.g. genesymbol=in.("CYP2D6","TPMT")).
"""
import json
import os
import time
import urllib.request

API = "https://api.cpicpgx.org/v1/"
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "cpic_api_raw")
os.makedirs(RAW, exist_ok=True)

GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "CYP2B6", "CYP3A5", "CYP4F2", "VKORC1", "DPYD",
         "TPMT", "NUDT15", "SLCO1B1", "UGT1A1", "G6PD", "HLA-A", "HLA-B", "IFNL3",
         "CFTR", "MT-RNR1", "RYR1", "CACNA1S"]
GENE_FILTER = "in.(" + ",".join('"%s"' % x for x in GENES) + ")"


def api(path, **params):
    url = API + path + ("?" + "&".join("%s=%s" % kv for kv in params.items()) if params else "")
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read())
        except Exception:
            if attempt == 2:
                raise
            time.sleep(3)


def save(name, data):
    with open(os.path.join(RAW, name + ".json"), "w") as fh:
        json.dump(data, fh, indent=1)
    print("%-24s %6d rows" % (name, len(data)))


def main():
    save("gene", api("gene", limit="5000"))
    save("pair", api("pair", limit="5000"))
    save("drug", api("drug", limit="5000"))
    save("guideline", api("guideline", limit="500"))
    save("recommendation", api("recommendation", limit="20000"))
    save("diplotype", api("diplotype", genesymbol=GENE_FILTER, limit="200000"))
    save("allele", api("allele", genesymbol=GENE_FILTER, limit="20000"))
    save("allele_definition", api("allele_definition", genesymbol=GENE_FILTER, limit="20000"))
    seqloc = api("sequence_location", genesymbol=GENE_FILTER, limit="200000")
    save("sequence_location", seqloc)

    ids = [s["id"] for s in seqloc]
    alv = []
    for i in range(0, len(ids), 150):
        alv += api("allele_location_value",
                   locationid="in.(" + ",".join(map(str, ids[i:i + 150])) + ")",
                   limit="200000")
    save("allele_location_value", alv)


if __name__ == "__main__":
    main()
