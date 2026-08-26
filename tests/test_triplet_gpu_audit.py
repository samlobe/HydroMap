from hydromap.engines.triplets import run_triplets_gpu


def test_tail_window_uses_explicit_physical_frame_interval() -> None:
    assert list(run_triplets_gpu._frame_indices(600, dt_ps=10.0, tail_ns=1.0, skip=1)) == list(
        range(500, 600)
    )


def test_sub_frame_tail_still_analyzes_final_frame() -> None:
    assert list(run_triplets_gpu._frame_indices(5, dt_ps=10.0, tail_ns=0.001, skip=1)) == [4]
