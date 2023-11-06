import click
import pickle
import pandas as pd

from glob import glob
from pathlib import Path
from scipy.stats import ks_2samp
from tart.utils.plotting import CoveragePlot
from tart.utils.read_parsing import AlignDat
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.activity_inference import Candidate
from tart.utils.filter_functions import default_check
from tart.utils.activity_inference import has_candidate_peak


def _log_cand_charac(align_charac: dict, cand: Candidate):
    align_charac["coverage_delta"] = cand.abs_cov_delta
    align_charac["coverage_delta_pval"] = cand.coverage_drop_pvalue
    align_charac["peak_summit"] = cand.summit
    align_charac["peak_width"] = cand.width
    # align_charac["peak_l/r"] = (cand.left, cand.right)
    align_charac["coverage_delta_relative"] = cand.rel_cov_delta
    align_charac["coverage_delta_noiseset"] = str(cand.coverage_delta_noise)

    align_charac["from_riboswith_end"] = cand.from_switch_end
    align_charac["from_riboswith_end_relative"] = cand.from_switch_end_relative

    align_charac["ks_stat"] = cand.symks_stat
    align_charac["ks_p"] = cand.symks_pval


def _process_candidate_list(
    candlist: list[Candidate],
    align_charac: dict,
    charac_local: list,
    ref: str,
    transcriptome: str,
):
    if len(candlist) == 0:
        align_charac["decision_note"] = "No suitable candidates"
        align_charac["decision"] = "fail"
        charac_local.append(align_charac)
        return "fail"

    for cand in candlist:
        align_charac = {}
        align_charac["rowid"] = ref
        align_charac["transcriptome"] = str(transcriptome)
        _log_cand_charac(align_charac, cand)
        decision = default_check(cand)
        align_charac["decision"] = decision

        if decision == "pass":
            charac_local.append(align_charac)
            return decision

    # Getting to this point means there were no passes
    align_charac = {}
    align_charac["rowid"] = ref
    align_charac["transcriptome"] = str(transcriptome)
    _log_cand_charac(align_charac, candlist[0])
    decision = default_check(candlist[0])
    align_charac["decision"] = decision
    charac_local.append(align_charac)

    return "fail"


@click.command()
@click.option(
    "-i",
    "--pick-root",
    required=True,
    help="Directory root for pickled outputs. Same directory as the parser output.",
)
@click.option(
    "-o",
    "--out-dir",
    required=True,
    help="Directory root for subdirectories with rendered filtered plots.",
)
@click.option(
    "--bin-size",
    default=10,
    show_default=True,
    help="Bin size for fragment ends binning.",
)
@click.option(
    "--min-cov-depth",
    default=50,
    show_default=True,
    help="Maximum read coverage (not inferred or clipped coverage) within the riboswitch region must be equal or greater than the specified threshold to proceed with pass/fail classification and include the alignment in the pass rate calculations.",
)
@click.option(
    "--ext-prop",
    nargs=2,
    default=(1.0, 1.0),
    show_default=True,
    help="Riboswitch size proportion outside the riboswitch region to include in the space searched for candidate peaks. \
        \
        For a 100bp riboswitch, passing 0.2 0.6 sets the search space from 20bp preceeding the riboswitch 5' end to 60bp beyond the riboswitch 3' end. \
        Similarly, passing -0.2 -0.6 sets the search space as 20bp into the 5' end and 60bp from the 3' end.",
)
@click.option(
    "--run-depr", is_flag=True, help="(Dev use) Run the deprecated version instead."
)
@click.option(
    "--conv",
    is_flag=True,
    help="(Dev use) Output plots with the ends convolution panel.",
)
def exec_main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop, run_depr, conv):
    if run_depr:
        raise ValueError("No deprecated function to run")

    else:
        main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop, conv)


def main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop, conv):
    # Determine MPI context
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

    # List of all pickled alignment data
    if rank == 0:
        total_files = glob(f"{pick_root}/**/*.p", recursive=True)

    else:
        total_files = None

    total_files = comm.bcast(total_files, root=0)

    mp_con.set_full_list(total_files)
    worker_list = mp_con.generate_worker_list()

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keep track of characteristics for each alignment
    charac_local = []

    # First pass; find suitable candidates if present
    for pick_file in worker_list:
        # Load in the alignment data object
        with open(pick_file, "rb") as f:
            alignDat: AlignDat = pickle.load(f)

        switch_size = alignDat.switch_end - alignDat.switch_start

        # Save path formation and sample/ref IDs
        ref = Path(pick_file[:-2]).name
        targetname = Path(pick_file).parts[-3]

        transcriptome = Path(pick_file).parts[-2]

        # Reset to confirm no carry-over between alignments ---------------
        align_charac = {}
        align_charac["rowid"] = ref
        align_charac["transcriptome"] = str(transcriptome)

        # Check whether there are enough reads to proceed.
        # This needs to be done to avoid clutter in the results
        # that failed the filter
        if (
            max(alignDat.readcov[alignDat.switch_start : alignDat.switch_end])
            < min_cov_depth
        ):
            align_charac["decision_note"] = "Minimum coverage not reached"
            charac_local.append(align_charac)
            continue

        lmarg, rmarg = ext_prop
        candlist: list[Candidate] = has_candidate_peak(
            alignDat, left_margin=lmarg, right_margin=rmarg
        )

        passfaildir = _process_candidate_list(
            candlist, align_charac, charac_local, ref, transcriptome
        )
        save_path = out_dir.joinpath(passfaildir, f"{ref}#{transcriptome}.png")
        save_path.parent.mkdir(exist_ok=True, parents=True)

        # Calculate and set info for binned raw ends
        alignDat.bin_rawends(bin_size=bin_size)

        # We still want to see a large region around the switch in the plots
        lplot = 1.0 if lmarg < 1.0 else lmarg
        rplot = 1.0 if rmarg < 1.0 else rmarg

        lbuff = int(switch_size * lplot)
        rbuff = int(switch_size * rplot)
        end_buffers = [lbuff, rbuff]

        pObj = CoveragePlot(alignDat, end_buffers)
        if conv:
            pObj._with_conv(save_path)

        else:
            pObj.default(save_path)

    charac_local_arr = comm.gather(charac_local, root=0)

    if rank == 0:
        if charac_local_arr is None:
            raise TypeError("Gather failed")
        characteristics = []
        for instance_arr in charac_local_arr:
            characteristics.extend(instance_arr)

        # Make dataframe
        # df = pd.DataFrame({"target_name": classes, "pass_rate": rates})
        # df.to_csv(f"{out_dir}/pass_rates.csv", index=False)
        df = pd.DataFrame(characteristics)
        df.to_csv(f"{out_dir}/characteristics.csv", index=False)
