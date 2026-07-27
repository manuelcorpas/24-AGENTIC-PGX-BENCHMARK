# tools/

Third-party binaries used by the revision analyses. Nothing here is committed:
the artefacts are large and belong to their upstream projects, so they are
fetched on demand and pinned by checksum.

## PharmCAT 3.4.0

Independent CPIC implementation used by `code/63-pharmcat-comparator.py` to
validate the executed mapping (revision item N2, Reviewer 2 points 2 and 7).

```bash
python code/63-pharmcat-comparator.py --fetch   # prints the exact commands
```

Expected artefact: `pharmcat-3.4.0-all.jar`, sha256
`9317ef632bf6c9786ff0d9d455d4c9f6d2882ebd66ad7256b4ae958ddf454741`.
Requires a JRE 17 or later. See `real-genome-arm/SOFTWARE.md`.

This README exists so the directory survives a clean clone: the comparator
resolves paths against it, and a missing directory would fail on a reviewer's
checkout while passing on ours. That is precisely the defect Reviewer 2 point 6
identified, and `tests/test_repo_paths.py` now catches it.
