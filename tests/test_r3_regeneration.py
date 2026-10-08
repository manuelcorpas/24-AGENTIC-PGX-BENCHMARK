from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_regeneration_commands_are_local_and_complete(tmp_path):
    s=spec_from_file_location('regen',ROOT/'docs/R3-2026-09-16/SOURCE/01-regenerate.py');m=module_from_spec(s);s.loader.exec_module(m)
    before=list(tmp_path.iterdir())
    commands=m.commands('python-test')
    assert len(commands)==4
    assert [Path(c[1]).name for c in commands]==['98-third-revision.py','65-hierarchical-stats.py','73-crossed-mixed-model.py','88-supp-normalisation-by-gene.py']
    assert all(c[0]=='python-test' for c in commands)
    assert list(tmp_path.iterdir())==before
