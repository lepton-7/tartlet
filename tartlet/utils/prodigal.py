import click
import pandas as pd

from pathlib import Path
from Bio import SeqIO, SeqRecord
from subprocess import run, PIPE
from tart.utils.utils import print, rowid
from typing import Any, Generator, Literal, Optional
from tart.utils.mpi_context import BasicMPIContext


class ORF:
    def __init__(self, record: SeqRecord.SeqRecord) -> None:
        """A prodigal ORF record from the translation output file (.faa)

        Args:
            record (SeqRecord.SeqRecord): _description_
        """
        self.strand: int = int(record.description.split(" # ")[3])
        self.start: int = int(record.description.split(" # ")[1])
        self.stop: int = int(record.description.split(" # ")[2])

        self.orf_from = self.stop if self.strand < 0 else self.start
        self.orf_to = self.start if self.strand < 0 else self.stop

        self.contig = "_".join(record.id.split("_")[:-1])
        # self.seq = record.seq


class ProdigalOutput:
    def __init__(self, filepath: str | Path) -> None:
        """Parser class to ingest prodigal translation output files (.faa)

        Args:
            filepath (str | Path): Prodigal .faa output filepath.
        """
        self.path = Path(filepath)

        self._orfs: list[ORF] = []
        for record in SeqIO.parse(self.path, "fasta"):
            self._orfs.append(ORF(record))

    def __iter__(self) -> Generator[ORF, Any, Any]:
        for orf in self._orfs:
            yield orf

    def find_downstream_orf(
        self, position: int, strand: Literal[-1, 1], contig: str
    ) -> ORF | None:
        """Find and return the first downstream ORF for given coordinates.

        Args:
            position (int): Find ORF beyond this point.
            strand (int): Search strand (1 or -1).
            contig (str): Search contig.

        Returns:
            ORF: First downstream ORF.
        """
        if strand > 0:
            for orf in self._orfs:
                if (
                    orf.contig == contig
                    and orf.strand == strand
                    and orf.orf_from > position
                ):
                    return orf

        else:
            for orf in reversed(self._orfs):
                if (
                    orf.contig == contig
                    and orf.strand == strand
                    and orf.orf_from < position
                ):
                    return orf

        # No downstream ORF for given position
        return None


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
            f"Successfully completed prodigal for {input_file.name} on worker {rank}."
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
        print(f"Started {mp_con.size} workers.")

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


@click.command()
@click.option(
    "--ledger",
    required=True,
    help="Path to the ledger containing riboswitch information.",
)
@click.option(
    "-i", "--prodigal-dir", required=True, help="Directory with prodigal outputs."
)
def record_orf_locations(ledger, prodigal_dir):
    ledger_path = Path(ledger)
    prodigal_dir = Path(prodigal_dir)

    df = pd.read_csv(ledger_path)

    # MPI setup for unique genomes
    mp_con = BasicMPIContext(list(pd.unique(df["genome_accession"])))
    worker_list = mp_con.generate_worker_list()

    local_entries = {}
    for genome_name in worker_list:
        genome_path = prodigal_dir.joinpath(f"{genome_name}.faa")

        if not genome_path.is_file:
            print(
                f"Prodigal output {genome_path} does not exist.\nSkipping for {genome_name}."
            )
            continue

        orfs = ProdigalOutput(genome_path)

        for i, row in df.iterrows():
            if row["genome_accession"] != genome_name:
                continue

            strand = -1 if row["strand"] == "-" else 1
            downstream = orfs.find_downstream_orf(
                row["seq_from"], strand, row["query_name"]
            )

            if downstream is not None:
                local_entries[rowid(row)] = [i, downstream.orf_from, downstream.orf_to]

    if mp_con.rank == 0:
        print("Gathering entries...")
    entries_arr = mp_con.comm.gather(local_entries, root=0)

    if mp_con.rank == 0 and entries_arr is not None:
        print("Completed gather.")
        # Check if ORF locations are already recorded and handle accordingly
        if "orf_from" in df.columns and "orf_to" in df.columns:
            print("ORF location columns already exist, overwriting empty entries.")
            from_list: list[str] = list(df["orf_from"])
            to_list: list[str] = list(df["orf_to"])

        else:
            from_list: list[str] = []
            to_list: list[str] = []
            for i in range(len(df)):
                from_list.append("")
                to_list.append("")

        # Consolidate
        entries = {}
        for inst_dict in entries_arr:
            entries.update(inst_dict)

        # Add entries to the ledger
        written = 0
        for k, v in entries.items():
            idx = v[0]
            f = v[1]
            t = v[2]

            if k != rowid(df.iloc[idx]):
                raise ValueError(
                    f"Row IDs did not match: {k} =/= {rowid(df.iloc[idx])}"
                )

            if not len(str(from_list[idx])) and not len(str(to_list[idx])):
                from_list[idx] = f
                to_list[idx] = t
                written += 1

        print(f"Recorded {written} downstream ORF locations.")
        df["orf_from"] = from_list
        df["orf_to"] = to_list

        df.to_csv(ledger_path, index=False)

    raise SystemExit(0)
