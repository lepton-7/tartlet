import click

from pathlib import Path
from typing import Optional
from subprocess import run, PIPE
from tart.utils.utils import print
from tart.utils.mpi_context import BasicMPIContext


def prodigal(
    input_file: Path,
    out_dir: Path,
    output_file: Optional[Path],
    trans_file: Optional[Path],
    options: tuple | list = [],
    rank: int = 0,
):
    # Coerce arguments to proper types:
    input_file = Path(input_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The command list passed to run() must contain only strings.
    options = [str(x) for x in options]

    input_file_name = input_file.stem

    # Set or coerce
    output_file = (
        out_dir.joinpath(f"{input_file_name}.gbk")
        if output_file is None
        else Path(output_file)
    )

    trans_file = (
        out_dir.joinpath(f"{input_file_name}.faa")
        if trans_file is None
        else Path(trans_file)
    )

    cmd = [
        "prodigal",
        "-i",
        f"{input_file}",
        "-o",
        f"{output_file}",
        "-a",
        f"{trans_file}",
        *options,
    ]
    call = run(cmd, stdout=PIPE, stderr=PIPE)

    if call.returncode:
        print(f"Prodigal failed on worker {rank}:")
        print(f"Command: {cmd}")
        print(f"Output: {call.stdout.decode('utf-8')}")
        print(f"Error: {call.stderr.decode('utf-8')}")

    else:
        print(
            f"Successfully completed prodigal for {input_file_name} on worker {rank}."
        )


@click.command()
@click.option(
    "-o", "--out-dir", required=True, help="Output directory for prodigal outputs."
)
@click.argument("total_files", nargs=-1)
def default_prodigal(out_dir, total_files: tuple or list):
    """Runs prodigal on the input files.

    Supports MPI acceleration.

    Args:\n
        out_dir (str): Directory for prodigal output files. Made if does not exist.
        total_files (tupleorlist): List or tuple of file paths to pass as prodigal inputs.
    """
    # MPI setup
    mp_con = BasicMPIContext([*total_files])
    worker_list = mp_con.generate_worker_list()

    if mp_con.rank == 0:
        print(f"Started {mp_con.size} workers")

    options = ["-q"]

    for fasta_path in worker_list:
        prodigal(
            input_file=fasta_path,
            out_dir=out_dir,
            options=options,
            rank=mp_con.rank,
            trans_file=None,
            output_file=None,
        )
