import click
import requests
import pandas as pd

from glob import glob
from pathlib import Path
from typing import Optional
from subprocess import run, PIPE
from tartlet.utils.mpi_context import BasicMPIContext
from tartlet.utils.utils import print, get_datapath_obj

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
    seq_file: str | Path,
    out_dir: str | Path,
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
@click.option("--is-list", is_flag=True, help="Argument is a path to a list")
@click.argument("total_files", nargs=-1)
def default_scan_for_riboswitches(
    out_dir, total_files: tuple | list, no_stats: bool, is_list: bool
):
    """Runs input files against the latest (14.9) rfam riboswitch covariance models.

    Supports MPI acceleration.

    Args:\n
        out_dir (str): Output directory for cmscan output files. Made if does not exist.
        total_files (tupleorlist): List or tuple of file paths to pass as cmscan inputs.
        no_stats (bool): Suppresses cmscan output.
    """
    # MPI setup

    if is_list:
        list_p = Path(total_files[0])
        print(f"Reading genome paths from {list_p}")
        with open(list_p, "r") as f:
            total_files = [x.rstrip() for x in f]

        print(f"Loaded in {len(total_files)} genomes")
    else:
        # Test if arg passed is a directory
        p = Path(total_files[0])
        if p.is_dir():
            print(f"Argument passed is a directory; looking for .fna in {p}.")
            total_files = glob(f"{p}/*.fna")

    mp_con = BasicMPIContext([*total_files])
    worker_list = mp_con.generate_worker_list()

    # Check if out_dir exists
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    if mp_con.rank == 0:
        print(f"Started {mp_con.size} workers.")

    for fasta_path in worker_list:
        riboswitch_cmscan(
            seq_file=fasta_path, out_dir=out_dir, rank=mp_con.rank, no_stats=no_stats
        )


def __infernal_to_df(fpath: Path, name: str = None):
    if name is None:
        name = fpath.stem.split(".")[0]
        print(f"Dataset name not provided for {fpath}, using {name}")

    rows = []
    # This is so that extra data that is not in the infernal table can be
    # added to the row data and parsed at the same time
    additionalCols = r" %s %s"

    with open(fpath, "r") as f:
        gen_acc = " "

        # infernal puts the query file (which has its genome accession ID) at
        # the end of the riboswitch tables, so the easy way to associate
        # genome IDs with their riboswitches is to just reverse the results
        # and append it to the results until the next table and a new genome ID is detected
        for line in reversed(f.readlines()):
            li = line.strip()

            # ignore everything except the results rows and add taxonomy data
            if not li.startswith("#"):
                rows.append(line.rstrip() + additionalCols % (gen_acc, name))

            # start of results from a new query genome
            if li.startswith("# Query file:"):
                # Extract the genome ID from the big query file path. Example:
                # Query file: ~/Emerge_MAGs_v1/Derep95/20110600_P1M_23.fna
                # first at every /, then take the last bit: 20110600_P1M_23.fna
                # and slice it to remove the ".fna".
                gen_acc = line.rstrip().split("/")[-1][:-4]

    table = []

    for row in rows:
        # Taking the strings from the infernal output and putting them into
        # a table to turn into a pd df.
        row_contents = row.split()

        # Because data in the "description of target" column has spaces,
        # the split() commands puts each word as its own element
        # into row_contents. To pull the relevant info out of the array,
        # slice and include in a way that doesnt touch
        # indices occupied by the description.
        table.append(
            [*row_contents[:16], *row_contents[-len(additionalCols.split(" ")) + 1 :]]
        )

    # Took most of this from the infernal results table header to pass onto the df
    col_lab = [
        "target_name",
        "accession",
        "query_name",
        "accession",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "e-value",
        "genome_accession",
        "dataset",
    ]

    df = pd.DataFrame(
        table,
        columns=col_lab,
    )

    # I dont think these columns are informative.
    df = df.drop(columns=["mdl_from", "mdl_to", "mdl"], axis="columns")

    return df

from os import getcwd

@click.command()
@click.option("-i", "--in-dir", required=True, help="Infernal results directory path")
@click.option("-o", "--out-path", required=True, help="Output .csv table path")
def parse_results(in_dir, out_path):
    in_dir = Path(in_dir)
    if not in_dir.is_dir():
        raise ValueError(f"Path {in_dir} is not a directory.")

    in_list = in_dir.glob(f"*.txt")
    print(f"Searched for files in: {in_dir}")

    dfs = []
    for fpath in in_list:
        print(f"Processing {fpath}")
        dfs.append(__infernal_to_df(fpath))

    if len(dfs) == 0:
        raise ValueError(f"No results files parsed. Is the input correct: {in_dir}?")

    cumulative = pd.concat(dfs, ignore_index=True)
    cumulative.to_csv(out_path, index=False)
