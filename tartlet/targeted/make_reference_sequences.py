# Goes through genomes in the given dset and records riboswitch nucleotide
# sequences to files categorised by riboswitch target_name

import os
import click
import pandas as pd
import multiprocessing as mp

from pathlib import Path
from Bio import SeqIO, Seq
from functools import partial
from collections import defaultdict
from tartlet.utils.mpi_context import BasicMPIContext
from tartlet.utils.utils import timestamped_print as print


def _get_buffered_seq(
    MAG_path: Path, table: pd.DataFrame, pre_delta: int, post_delta: int
) -> list[tuple[str, dict]]:
    """Called in an MP pool"""

    MAGDict = {x.id: str(x.seq) for x in SeqIO.parse(MAG_path, "fasta")}
    subset = table[table["genome_accession"] == MAG_path.stem]
    ret = []

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
        classname: str = row["target_name"]
        qname = row["query_name"]
        frm = row["seq_from"]
        to = row["seq_to"]
        strand = row["strand"]

        # No spaces to ensure the entire string is recognised as the ID
        rowid = f"{classname}#{qname}#{frm}#{to}#{strand}"

        ret.append((classname, {rowid: str(switch)}))
        # seqs_local[classname].update({rowid: switch})
    return ret


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
    help="Dump alignment reference sequences in fasta format to this directory. Will be created if it does not exist.",
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
    help="Number of nucleotides upstream of riboswitches to capture in the generated reference sequence.",
)
@click.option(
    "--post-del",
    "post_delta",
    default=500,
    show_default=True,
    help="Number of nucleotides downstream of riboswitches to capture in the generated reference sequence.",
)
@click.option(
    "--unify",
    is_flag=True,
    help="(Dev use) Write all sequences to single file.",
)
def main(ledger_path, out_dir, genome_dir, dset, pre_delta, post_delta, unify):
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

    if rank == 0:
        print(f"Started {size} worker(s) with {alloc_cpus} cores each")

    local_path_list = mp_con.generate_worker_list()

    # Setup the local dictionary that stores the sequences

    # Should have structure: {
    #   "SAM" : {"SAM # 3300037155_10_Ga0395911_000841 # 6843 # 6736 # -" : "ATGC...", ...
    #   },
    #   "TPP" : {" TPP # ...": "GCTAGATT...", ...
    #   }, ...
    # }
    seqs_local = defaultdict(dict)
    # seqs_local = {str(rank) : rank * 50}

    if mp_con.is_active:
        with mp.Pool(alloc_cpus) as pool:
            print(f"Started a pool with {alloc_cpus} processes") if rank == 0 else None
            subset_f = partial(
                _get_buffered_seq,
                table=table,
                pre_delta=pre_delta,
                post_delta=post_delta,
            )
            ret_l = pool.map(subset_f, local_path_list)

        # Collect all the buffered sequences into the local dictionary
        _ = [seqs_local[classname].update(d) for tup in ret_l for classname, d in tup]

    seqs_arr = comm.gather(seqs_local, root=0)

    # defaultdict makes merging worker dicts trivial
    seqs_ledger = defaultdict(dict)

    # Consolidate on root thread
    if rank == 0:
        if seqs_arr is None:
            raise ValueError("Gather failed.")

        print(f"Completed gather on 0")

        for instance_dict in seqs_arr:
            for switchclass, seqs in instance_dict.items():
                seqs_ledger[switchclass].update(seqs)

        # Make the sub directory to save sequences in fasta format
        out_dir.mkdir(parents=True, exist_ok=True)

        num_switch_classes = len(seqs_ledger)
    else:
        num_switch_classes = None

    # broadcast the number of distinct switch class dicts in the ledger
    # so the extra threads can exit
    if not unify:
        num_switch_classes = comm.bcast(num_switch_classes, root=0)
    else:
        num_switch_classes = 0

    def multithreaded_writeout(num_switch_classes, seqs_ledger):
        # Setup a dummy list to be able to subscript in the for loop and
        # send each switch class to one worker thread
        if not rank:
            dummy_ledger = [(key, val) for key, val in seqs_ledger.items()]

        # Iterate over sub dictionaries
        for idx in range(num_switch_classes):
            # Should probably refactor to use scatter instead of individual sends
            if not rank:
                key, val = dummy_ledger[idx]  # type: ignore

                comm.send(key, dest=idx + 1, tag=10)  # type: ignore
                comm.send(val, dest=idx + 1, tag=100)  # type: ignore

            elif rank == idx + 1:
                classname = comm.recv(source=0, tag=10)  # type: ignore
                sub_d = comm.recv(source=0, tag=100)  # type: ignore

                # Write riboswitch sequences to disk
                write_step(classname, sub_d, out_dir)  # type: ignore

    # Each riboswitch class sub-dictionary is sent to a worker for disk writes
    # if there are enough workers for that. Otherwise the root performs the write out

    if unify:
        unified_writeout(seqs_ledger, out_dir) if rank == 0 else None

    elif not unify and size > int(num_switch_classes):
        # Exit workers that are not needed
        if rank > int(num_switch_classes):
            raise SystemExit(0)

        if rank > 0:
            classname = None
            sub_d = None

        # TODO: Refactor to use both MPI and mp
        multithreaded_writeout(num_switch_classes, seqs_ledger)

    else:
        singlethreaded_writeout(seqs_ledger, out_dir) if rank == 0 else None


def write_step(classname: str, sub_d: dict, out_dir: Path):
    # Write riboswitch sequences to disk
    # if rank > 0:
    fpath = out_dir.joinpath(f"{classname}.fna")

    with open(fpath, "w") as f:
        for key, val in sub_d.items():
            f.write(f">{key}\n")
            f.write(f"{val}\n")


def singlethreaded_writeout(seqs_ledger: dict, out_dir: Path):
    for classname, sub_d in seqs_ledger.items():
        write_step(classname, sub_d, out_dir)


def unified_writeout(seqs_ledger: dict, out_dir: Path):
    print("Starting unified write-out.")
    unified_dict = {}
    for _, sub_d in seqs_ledger.items():
        unified_dict.update(sub_d)

    write_step("unified", unified_dict, out_dir)
