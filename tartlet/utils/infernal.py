from typing import Optional
import click
import requests

from pathlib import Path
from subprocess import run, PIPE
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.utils import print, get_datapath_obj

current_ver = 14.9

# As of rfam 14.9
rfam_riboswitch_accessions = [
    "RF00050",
    "RF00059",
    "RF00080",
    "RF00162",
    "RF00167",
    "RF00168",
    "RF00174",
    "RF00234",
    "RF00379",
    "RF00380",
    "RF00442",
    "RF00504",
    "RF00521",
    "RF00522",
    "RF00634",
    "RF01054",
    "RF01055",
    "RF01056",
    "RF01057",
    "RF01482",
    "RF01510",
    "RF01689",
    "RF01704",
    "RF01725",
    "RF01727",
    "RF01734",
    "RF01739",
    "RF01750",
    "RF01767",
    "RF01786",
    "RF01826",
    "RF01831",
    "RF02680",
    "RF02683",
    "RF02885",
    "RF02912",
    "RF03057",
    "RF03058",
    "RF03071",
    "RF03072",
]


def make_riboswitch_cm(ver: float = current_ver):
    # Get a pathlib.Path object for the data subdirectory
    data_dir_obj = get_datapath_obj()
    cm = str(data_dir_obj.joinpath(f"riboswitches_{ver}.cm"))

    with open(cm, "w") as f:
        for acc in rfam_riboswitch_accessions:
            address = f"https://rfam.org/family/{acc}/cm"

            r = requests.get(address)

            f.write(r.text)
            print(f"Ingested {acc}")

    print("Starting cmpress")

    run(["cmpress", "-F", f"{cm}"])


# DEPRECATE: not necessary anymore
def make_clanin(ver: float = current_ver):
    url = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/{ver}/Rfam.clanin"

    # Get a pathlib.Path object for the data subdirectory
    data_dir_obj = get_datapath_obj()
    clanin = str(data_dir_obj.joinpath(f"rfam_{ver}.clanin"))

    with open(clanin, "w") as f:
        r = requests.get(url)
        f.write(r.text)


def cmscan(
    seq_file: Path,
    out_dir: Path,
    cm_path: Path,
    options: tuple | list,
    out_file: Optional[Path] = None,
    rank: int = 0,
    no_stats: bool = False,
):
    # Coerce arguments to proper types:
    seq_file = Path(seq_file)
    out_dir = Path(out_dir)
    cm_path = Path(cm_path)

    # The command list passed to run() must contain only strings.
    options = [str(x) for x in options]

    seq_file_name = seq_file.name

    # Set or coerce
    out_file = (
        out_dir.joinpath(f"{seq_file_name}.txt") if out_file is None else Path(out_file)
    )

    cmd = [
        "cmscan",
        "--cut_ga",
        "--rfam",
        "--nohmmonly",
        "--noali",
        "--tblout",
        f"{out_file}",
        *options,
        f"{cm_path}",
        f"{seq_file}",
    ]
    call = run(cmd, stdout=PIPE, stderr=PIPE)

    if call.returncode:
        print(f"Failed on worker {rank}:")
        print(f"Command: {cmd}")
        print(f"Output: {call.stdout.decode('utf-8')}")
        print(f"Error: {call.stderr.decode('utf-8')}")

    else:
        if not no_stats:
            print(call.stdout.decode("utf-8"))
        print(
            f"Successfully completed cmscan for {seq_file_name} against {cm_path.name} on worker {rank}."
        )


def riboswitch_cmscan(
    seq_file: str or Path,
    out_dir: str or Path,
    options: tuple | list = [],
    ver: float = current_ver,
    **kwargs,
):
    """Runs a cmscan instance against a riboswitch-specific covariance model.

    Args:
        seq_file (strorPath): cmscan input file.
        out_dir (strorPath): Directory for cmscan output.
        options (tupleorlist, optional): Optional options to pass to the cmscan call. Defaults to [].
        ver (float, optional): CM rfam version. Defaults to current_ver.

    **kwargs:
        Inserted into the cmscan() call.

    Raises:
        FileNotFoundError: Riboswitch CM not found.
    """
    # Coerce input
    ver = float(ver)

    data_dir_obj = get_datapath_obj()
    switch_cm = data_dir_obj.joinpath(f"riboswitches_{ver}.cm")

    if not switch_cm.exists():
        raise FileNotFoundError(
            f"Riboswitch CM for RFAM {ver} does not exist. Call make_riboswitch_cm({ver})."
        )

    cmscan(
        Path(seq_file),
        Path(out_dir),
        switch_cm,
        options=[*options],
        **kwargs,
    )


@click.command()
@click.option(
    "-o", "--out-dir", required=True, help="Output directory for cmscan output."
)
@click.option("--no-stats", is_flag=True, help="Supresses cmscan output to the console")
@click.argument("total_files", nargs=-1)
def default_scan_for_riboswitches(out_dir, total_files: tuple or list, no_stats: bool):
    """Runs input files against the latest (14.9) rfam riboswitch covariance models.

    Supports MPI acceleration.

    Args:\n
        out_dir (str): Output directory for cmscan output files.
        total_files (tupleorlist): List or tuple of file paths to pass as cmscan inputs.
        no_stats (bool): Suppresses cmscan output.
    """
    # MPI setup
    mp_con = BasicMPIContext([*total_files])
    worker_list = mp_con.generate_worker_list()

    print(f"Started {mp_con.size} workers")

    for fasta_path in worker_list:
        riboswitch_cmscan(
            seq_file=fasta_path, out_dir=out_dir, rank=mp_con.rank, no_stats=no_stats
        )
