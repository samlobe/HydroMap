from pathlib import Path

from hydromap.manifest import sha256_file, stable_payload_hash, summarize_results_csv


def test_sha256_file_is_stable(tmp_path: Path) -> None:
    f = tmp_path / "a.txt"
    f.write_text("hydromap\n", encoding="utf-8")

    h1 = sha256_file(f)
    h2 = sha256_file(f)
    assert h1 == h2


def test_results_summary(tmp_path: Path) -> None:
    csv_path = tmp_path / "results.csv"
    csv_path.write_text(
        "MDAnalysis_selection_strings,Fdewet_pred,U_pw\n"
        "resid 1,4.0,-120.0\n"
        "resid 2,6.0,-90.0\n",
        encoding="utf-8",
    )
    summary = summarize_results_csv(csv_path)
    assert summary["rows"] == 2
    assert summary["fdewet_mean"] == 5.0
    assert summary["fdewet_min"] == 4.0
    assert summary["fdewet_max"] == 6.0


def test_stable_payload_hash_order_independent() -> None:
    a = {"x": 1, "y": [1, 2, 3]}
    b = {"y": [1, 2, 3], "x": 1}
    assert stable_payload_hash(a) == stable_payload_hash(b)
