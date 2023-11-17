import click
import pickle
import tarfile as tf

from shutil import rmtree
from glob import glob
from pathlib import Path
from tart.utils.plotting import CoveragePlot
from tart.utils.read_parsing import SortedBAM
from tart.utils.mpi_context import BasicMPIContext


@click.command()
@click.option(
    "-i",
    "--bam-dir",
    required=True,
    help="Finds sorted BAM files within its subdirectories.",
)
@click.option(
    "-o",
    "--out-dir",
    required=True,
    help="Directory root to place parsed output subdirectories.",
)
@click.option(
    "--bounds-file",
    required=True,
    help="File that stores the riboswitch bounds map used during visualisation.",
)
@click.option(
    "--min-coverage",
    default=8,
    show_default=True,
    help="Minimum converage threshold (in reads). Coverage for at least one position across the reference must be above the value specified, otherwise alignment output is not generated.",
)
@click.option(
    "--plots",
    "outPlots",
    is_flag=True,
    help="Render outputs via matplotlib and save.",
)
@click.option(
    "--picks",
    "outPickles",
    is_flag=True,
    help="Save outputs as binary pickles.",
)
@click.option(
    "--allow-soft-clips",
    is_flag=True,
    help="Consider the start/end of the soft clipped regions when determining where the likely fragment ends are from reads. If not set, only consider the aligned portions of the reads regardless of any soft clipped regions extending beyond.",
)
def main(
    bam_dir, out_dir, bounds_file, min_coverage, outPlots, outPickles, allow_soft_clips
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # List of all sorted BAMS to process
    total_files = glob(f"{bam_dir}/**/*.sorted.bam")

    # MPI setup
    mp_con = BasicMPIContext(total_files)
    worker_list = mp_con.generate_worker_list()

    for bam_path in worker_list:
        bam_wrap = SortedBAM(bam_path, bounds_file)

        # This looks something like this:
        # /users/PDS0325/sachitk26/ribo_acclim/metatranscriptomics/
        # switch_alignment/switch_seqs_delta500/alignments_full/AdoCbl_riboswitch/AdoCbl_riboswitch.MainAutochamber.201707_E_2_20to24.sorted.bam
        bam_path = Path(bam_path)
        bam_split = bam_path.parts

        # Looks like:
        # save_dir/AdoCbl_riboswitch/AdoCbl_riboswitch.MainAutochamber.201707_E_2_20to24
        #
        # The [:-11] removes the .sorted.bam suffix from the path, but Path objects are
        # not subscriptable, so it needs to be typecasted to str, sliced, then converted
        # back to a Path. Don't @ me this works just fine.
        save_dir = Path(str(out_dir.joinpath(*bam_split[-2:]))[:-11])

        # Make sure the directory exists; create if not.
        save_dir.mkdir(parents=True, exist_ok=True)

        alignDat_arr = bam_wrap.generate_ref_alignment_data(allow_soft_clips)

        for alignDat in alignDat_arr:
            if outPlots and alignDat.is_coverage_threshold("read", min_coverage):
                save_path = save_dir.joinpath(f"{alignDat.ref}.png")
                CoveragePlot(alignDat, [40, 40]).default(save_path)

            if outPickles and alignDat.is_coverage_threshold("read", min_coverage):
                save_path = save_dir.joinpath(f"{alignDat.ref}.p")
                with open(save_path, "wb") as f:
                    pickle.dump(alignDat, f)

    # This is just so that the root waits until all the workers are done
    done_workers = mp_con.comm.gather(mp_con.rank, root=0)

    if mp_con.rank == 0 and done_workers is not None:
        if len(done_workers) == mp_con.size:
            tarpath = out_dir.parent.joinpath(f"{out_dir.name}.tar.gz")
            with tf.open(tarpath, "w:gz") as picktar:
                print(f"Archiving pickled data")
                picktar.add(out_dir)

            print(f"Removing directory")
            rmtree(out_dir)

        else:
            print(f"Did all workers not finish?: {done_workers}")

    raise SystemExit(0)
