from __future__ import annotations

import os


def compute_cpu_workers(max_cpu_workers: int, reserve_cpus: int, profile: str) -> int:
    total = os.cpu_count() or 1
    budget = total - reserve_cpus
    if budget < 1:
        budget = 1

    hard_cap = max(1, min(max_cpu_workers, budget))

    if profile == "gpu_fast":
        return hard_cap

    if profile == "balanced":
        conservative = max(1, budget // 2)
        return max(1, min(hard_cap, conservative))

    return hard_cap
