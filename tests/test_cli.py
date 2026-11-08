from argparse import Namespace

from onco_dynamics_x.cli import run_simulation


def test_cli_can_run_small_3d_simulation(tmp_path):
    output = tmp_path / "summary.csv"
    args = Namespace(
        grid_size=7,
        dimensions=3,
        dt=0.0001,
        steps=2,
        save_every=1,
        output=str(output),
    )
    run_simulation(args)
    assert output.exists()
    assert len(output.read_text(encoding="utf-8").splitlines()) == 4
