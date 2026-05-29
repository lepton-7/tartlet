# Goes through derep95 MAGs, finds riboswitch sequences with a buffer around the sequence
# and returns a dict for indices of the full switch + buffer sequences in the contig

import os
import json
import click
import pandas as pd
import multiprocessing as mp

from pathlib import Path
from Bio import SeqIO, Seq
from functools import partial
from collections import defaultdict
from tartlet.utils.mpi_context import BasicMPIContext
from tartlet.utils.utils import timestamped_print as print


def _get_rowid_bounds(
    MAG_path: Path, table: pd.DataFrame, pre_delta: int, post_delta: int
) -> list[dict]:
    MAGDict = {x.id: str(x.seq) for x in SeqIO.parse(MAG_path, "fasta")}
    subset = table[table["genome_accession"] == MAG_path.stem]
    ret_l = []

    for _, row in subset.iterrows():  # iterate over MAG riboswitches
        # Infernal start and stop entries are relative to canonical 5' -> 3';
        # need to account for that when slicing the sequence
        if row["strand"] == "+":
            start = int(row["seq_from"])
            end = int(row["seq_to"])

        elif row["strand"] == "-":
            start = int(row["seq_to"])
            end = int(row["seq_from"])

        else:
            print(f"Yikes: strand notation not recognised.")  # should not be possible
            continue

        contigseq = MAGDict[row["query_name"]]

        # Adjust bounds with delta
        start -= pre_delta
        end += post_delta

        # Validate bounds
        if start < 0:
            start = 0
        if end > len(contigseq):
            end = len(contigseq)

        # Find query seq
        switch = contigseq[start:end]

        # switch is on the (-) strand
        if row["strand"] == "-":
            switch = Seq.reverse_complement(switch)

        # Add the sequence to the dict
        classname = row["target_name"]
        qname = row["query_name"]
        frm = row["seq_from"]
        to = row["seq_to"]
        strand = row["strand"]

        o_from = row["orf_from"]
        o_to = row["orf_to"]

        # No spaces to ensure the entire string is recognised as the ID
        rowid = f"{classname}#{qname}#{frm}#{to}#{strand}"

        ret_l.append({rowid: {"ref": (start, end), "orf": (o_from, o_to)}})

    return ret_l


@click.command()
@click.option(
    "--ledger",
    "ledger_path",
    required=True,
    help="Path to the ledger containing riboswitch information.",
)
@click.option(
    "--out-dir",
    required=True,
    help="Dump json map to this directory. Will be created if it does not exist.",
)
@click.option(
    "--genome",
    "genome_dir",
    required=True,
    help="Path to the genome/genome(s) directory. Genomes must have the .fna prefix. If path is to a directory, all genomes within will be processed together unless ledger rows are selected using --dset.",
)
@click.option(
    "--dset",
    default=None,
    show_default=True,
    help="Select specific rows in the ledger according to the dataset column. Handy if there are multiple genomes in --genome-dir and only a subset need to have reference sequences generated in a group.",
)
@click.option(
    "--pre-del",
    "pre_delta",
    default=500,
    show_default=True,
    help="Number of nucleotides upstream of riboswitches captured in the generated reference sequence.",
)
@click.option(
    "--post-del",
    "post_delta",
    default=500,
    show_default=True,
    help="Number of nucleotides downstream of riboswitches captured in the generated reference sequence.",
)
def main(ledger_path, out_dir, genome_dir, dset, pre_delta, post_delta):
    # Read in the results table
    table = pd.read_csv(ledger_path)

    # Drop all but the relevant columns
    table = table[
        [
            "target_name",
            "query_name",
            "seq_from",
            "seq_to",
            "strand",
            "trunc",
            "genome_accession",
            "dataset",
            "orf_from",
            "orf_to",
        ]
    ]
    if dset is not None:
        table = table[table["dataset"] == dset]

    genome_dir = Path(genome_dir)
    out_dir = Path(out_dir)

    if not genome_dir.exists():
        raise ValueError(f" File/directory {genome_dir} does not exist")

    if genome_dir.is_dir():
        genomes_list = [*genome_dir.glob(f"*.fna")]

    elif genome_dir.is_file():
        genomes_list = [genome_dir]

    else:
        raise ValueError(
            f"Genome directory '{genome_dir} is neither a directory nor a file'"
        )

    # Setup MPI to parse switch sequences
    mp_con = BasicMPIContext(genomes_list)
    comm = mp_con.comm
    size = mp_con.size
    rank = mp_con.rank

    alloc_cpus = len(os.sched_getaffinity(0))

    (
        print(f"Started {size} worker(s) with {alloc_cpus} cores each")
        if rank == 0
        else None
    )

    local_path_list = mp_con.generate_worker_list()

    # Setup the local dictionary that stores the sequences
    seqs_local = defaultdict()
    if mp_con.is_active:
        with mp.Pool(alloc_cpus) as pool:
            print(f"Started a pool with {alloc_cpus} processes") if rank == 0 else None
            subset_f = partial(
                _get_rowid_bounds,
                table=table,
                pre_delta=pre_delta,
                post_delta=post_delta,
            )
            ret_l = pool.map(subset_f, local_path_list)
        _ = [seqs_local.update(d) for subl in ret_l for d in subl]

    seqs_arr = comm.gather(seqs_local, root=0)

    # Consolidate on root thread
    if rank == 0:
        if seqs_arr is None:
            raise ValueError("Gather failed.")
        print("Completed gather on 0")
        # defaultdict makes merging worker dicts trivial
        seqs_ledger = defaultdict()

        for instance_dict in seqs_arr:
            seqs_ledger.update(instance_dict)

        out_dir.mkdir(parents=True, exist_ok=True)

        with open(out_dir.joinpath(f"rowid_to_bounds.json"), "w") as f:
            json.dump(seqs_ledger, f, indent=1)

    else:
        raise SystemExit(0)
