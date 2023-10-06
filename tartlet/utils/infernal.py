import requests
import click

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
    options: tuple or list,
    out_file: Path = None,
):
    seq_file_name = seq_file.name

    if out_file is None:
        out_file = out_dir.joinpath(f"{seq_file_name}.txt")

    call = run(
        [
            "cmscan",
            "--cut_ga",
            "--rfam",
            "--nohmmonly",
            "--noali",
            "--tblout",
            f"{out_file}",
            f"{cm_path}",
            f"{seq_file}",
            *options,
        ],
        stdout=PIPE,
        stderr=PIPE,
    )
    if call.returncode:
        print("Failed:\n", call.stderr.decode("utf-8"))
    else:
        print(call.stdout.decode("utf-8"))
        print(
            f"Successfully completed cmscan on {seq_file_name} against {cm_path.name}."
        )


def riboswitch_cmscan(
    seq_file: str or Path,
    out_dir: str or Path,
    options: tuple or list = [],
    ver: float = current_ver,
):
    data_dir_obj = get_datapath_obj()
    switch_cm = data_dir_obj.joinpath(f"riboswitches_{ver}.cm")
    clanin = data_dir_obj.joinpath(f"rfam_{ver}.clanin")

    if not switch_cm.exists():
        raise ValueError(
            f"Riboswitch CM for RFAM {ver} does not exist. Call make_riboswitch_cm({ver})."
        )

    if not clanin.exists():
        raise ValueError(
            f"Clanin file for RFAM {ver} does not exist. Call make_clanin({ver})."
        )

    cmscan(
        Path(seq_file),
        Path(out_dir),
        switch_cm,
        options=["--clanin", clanin, *options],
    )


@click.command()
@click.option(
    "-o", "--out-dir", required=True, help="Output directory for cmscan output."
)
@click.argument("total_files", nargs=-1)
def dafault_scan_for_riboswitches(out_dir, total_files: tuple or list):
    """Runs input files against the latest (14.9) rfam riboswitch covariance models.

    Supports MPI acceleration.

    Args:
        out_dir (str): Output directory for cmscan output files.
        total_files (tupleorlist): List or tuple of file paths to pass as cmscan inputs.
    """
    # MPI setup
    mp_con = BasicMPIContext([*total_files])
    worker_list = mp_con.generate_worker_list()

    for fasta_path in worker_list:
        riboswitch_cmscan(Path(fasta_path), Path(out_dir))
