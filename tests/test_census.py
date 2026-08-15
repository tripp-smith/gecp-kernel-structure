from pathlib import Path

from kernelgecp.census import run_census


def test_census_is_byte_reproducible(tmp_path: Path) -> None:
    config = tmp_path / "census.toml"
    config.write_text(
        """schema_version = 1
cutoffs = [1, 4]
tolerances = [1e-6]
precision_bits = 128
max_rank = 12
grid_order = 6
pivot = "grid"
""",
        encoding="utf-8",
    )
    root = Path(__file__).parents[1]
    first = tmp_path / "first.jsonl"
    second = tmp_path / "second.jsonl"
    run_census(config, first, root=root)
    run_census(config, second, root=root)
    assert first.read_bytes() == second.read_bytes()
    assert len(first.read_text(encoding="utf-8").splitlines()) == 2
