from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_audit_contains_all_candidates_and_explicit_scoring_limits():
    s=spec_from_file_location('audit',ROOT/'docs/R3-2026-09-16/SOURCE/02-build-audit.py');m=module_from_spec(s);s.loader.exec_module(m)
    text=m.candidate_audit()
    assert text.count('## Case: ')==8
    assert '95.73%' not in text
    assert 'not independent clinical adjudication' in text
    assert '0/24 (0.00)' in text
    assert 'https://doi.org/10.1002/cpt.1319' in text
    assert text.count('```text')==16
