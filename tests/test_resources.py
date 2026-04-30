from hydromap.resources import compute_cpu_workers


def test_gpu_fast_uses_hard_cap(monkeypatch) -> None:
    monkeypatch.setattr("os.cpu_count", lambda: 16)
    assert compute_cpu_workers(max_cpu_workers=8, reserve_cpus=2, profile="gpu_fast") == 8


def test_balanced_is_more_conservative(monkeypatch) -> None:
    monkeypatch.setattr("os.cpu_count", lambda: 16)
    assert compute_cpu_workers(max_cpu_workers=16, reserve_cpus=2, profile="balanced") == 7


def test_compute_cpu_workers_never_below_one(monkeypatch) -> None:
    monkeypatch.setattr("os.cpu_count", lambda: 2)
    assert compute_cpu_workers(max_cpu_workers=8, reserve_cpus=5, profile="balanced") == 1
