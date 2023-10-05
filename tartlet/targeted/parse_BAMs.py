import pysam
import click
import json, pickle

from glob import glob
from pathlib import Path
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.parsefuncs import generate_plot_data, plot_gen


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
    help="Minimum converage threshold (in reads). Coverage for at least one position across the reference must be above the value specified, otherwise output is not generated for this alignment.",
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
def main(bam_dir, out_dir, bounds_file, min_coverage, outPlots, outPickles):
    # bam_dir = argv[1]

    # out_dir = Path(argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    # Set these to determine outputs ----------------------------------------------
    # min_coverage = 8

    # outPlots = False
    # outPickles = True
    # -----------------------------------------------------------------------------

    # Import the dict that stores MAG contig positions of the alignment references
    with open(f"{bounds_file}", "r") as f:
        refBounds = json.load(f)

    # List of all sorted BAMS to process
    total_files = glob(f"{bam_dir}/**/*.sorted.bam")

    # MPI setup
    mp_con = BasicMPIContext(total_files)
    worker_list = mp_con.generate_worker_list()

    for bam_path in worker_list:
        bam = pysam.AlignmentFile(bam_path, "rb")

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

        outDict = generate_plot_data(bam, refBounds)

        for ref, alignTup in outDict.items():
            start, end = alignTup[2][0], alignTup[2][1]
            readcov = alignTup[0][0]

            if outPlots and max(readcov[start:end]) >= min_coverage:
                save_path = save_dir.joinpath(f"{ref}.png")
                plot_gen(ref, alignTup, save_path, buff=40)

            if outPickles and max(readcov[start:end]) >= min_coverage:
                save_path = save_dir.joinpath(f"{ref}.p")
                with open(save_path, "wb") as f:
                    pickle.dump(alignTup, f)
