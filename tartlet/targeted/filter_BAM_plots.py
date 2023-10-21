import click
import pickle
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.activity_inference import (
    plot_gen,
    is_interesting,
    has_candidate_peak,
)
from tart.utils.activity_inference import Candidate
from tart.utils.read_parsing import AlignDat


def _log_cand_charac(align_charac: dict, cand: Candidate):
    align_charac["coverage_delta"] = cand.abs_cov_delta
    align_charac["coverage_delta_relative"] = cand.rel_cov_delta
    align_charac["coverage_delta_noiseset"] = str(cand.coverage_delta_noise)

    align_charac["from_riboswith_end"] = cand.from_switch_end
    align_charac["from_riboswith_end_relative"] = cand.from_switch_end_relative


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
    "--run-depr", is_flag=True, help="(Dev use) Run the deprecated version instead."
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
def exec_main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop, run_depr):
    if run_depr:
        depr_main(pick_root, out_dir, bin_size)

    else:
        main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop)


def main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop):
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

        # Reset to confirm no carry-over between alignments
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
        cand: Candidate = has_candidate_peak(
            alignDat, left_margin=lmarg, right_margin=rmarg
        )
        passfaildir = "fail" if cand is None else "pass"
        align_charac["decision"] = passfaildir
        save_path = out_dir.joinpath(passfaildir, f"{transcriptome}#{ref}.png")

        if cand is None:
            align_charac["decision_note"] = "No suitable candidates"
            charac_local.append(align_charac)
            continue

        # Extract candidate characteristics and organise in-place
        _log_cand_charac(align_charac, cand)
        charac_local.append(align_charac)

        save_path.parent.mkdir(exist_ok=True, parents=True)

        # Calculate and set info for binned raw ends
        alignDat.bin_rawends(bin_size=bin_size)

        lbuff = int(switch_size * lmarg)
        rbuff = int(switch_size * rmarg)

        # plot_gen(
        #     ref,
        #     alignTup,
        #     str(save_path),
        #     bin_size=bin_size,
        #     bin_ax=bin_ax,
        #     lbuff=lbuff,
        #     rbuff=rbuff,
        # )

    charac_local_arr = comm.gather(charac_local, root=0)

    if rank == 0:
        characteristics = []
        for instance_arr in charac_local_arr:
            characteristics.extend(instance_arr)

        # Make dataframe
        # df = pd.DataFrame({"target_name": classes, "pass_rate": rates})
        # df.to_csv(f"{out_dir}/pass_rates.csv", index=False)
        df = pd.DataFrame(characteristics)
        df.to_csv(f"{out_dir}/characteristics.csv", index=False)


def depr_main(pick_root, out_dir, bin_size):
    # Determine MPI context
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

    # List of all pickled plots
    if rank == 0:
        total_files = glob(f"{pick_root}/**/*.p", recursive=True)

    else:
        total_files = None

    total_files = comm.bcast(total_files, root=0)

    mp_con.set_full_list(total_files)
    worker_list = mp_con.generate_worker_list()

    # kernel = gen_kernel(kernel_size=41, std_dev=5)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keeps track of which plots filtered through
    pass_rate_local = defaultdict(list)

    for pick_file in worker_list:
        with open(pick_file, "rb") as f:
            # ([read cov, inferred frag cov, clipped cov], ends, (swtch start, end))
            cov, ends, (switch_start, switch_end) = pickle.load(f)
            readcov, infercov, clipcov = cov

        # conv_ends = np.convolve(ends, kernel, "same")

        # convalignTup = (cov, conv_ends, (switch_start, switch_end))
        alignTup = (cov, ends, (switch_start, switch_end))

        # Check whether there are enough reads to proceed.
        # This needs to be done to avoid clutter in the results
        # that failed the filer
        if max(readcov) < 50:
            continue

        isActive = is_interesting(alignTup, windowfrac=0.2, threshtol=0.2)
        passfaildir = "pass" if isActive else "fail"

        # Save path formation
        ref = Path(pick_file[:-2]).name

        # Riboswitch class name
        targetname = Path(pick_file).parts[-3]

        # Tally pass vs fail
        pass_rate_local[ref].append(isActive)

        core_sample = Path(pick_file).parts[-2]
        save_path = out_dir.joinpath(passfaildir, f"{core_sample}#{ref}.png")

        save_path.parent.mkdir(exist_ok=True, parents=True)

        alignTup, bin_ax = bin_counts(alignTup, bin_size=bin_size)

        plot_gen(
            ref,
            alignTup,
            str(save_path),
            bin_size=bin_size,
            bin_ax=bin_ax,
        )

    charac_local_arr = comm.gather(pass_rate_local, root=0)

    if rank == 0:
        pass_rates = defaultdict(list)
        for instance_dict in charac_local_arr:
            for key, val in instance_dict.items():
                pass_rates[key].extend(val)

        # Output pass rates to file
        classes = []
        rates = []
        for key, val in pass_rates.items():
            classes.append(key)
            rates.append(sum(val) / len(val))

        df = pd.DataFrame({"target_name": classes, "pass_rate": rates})
        df.to_csv(f"{out_dir}/pass_rates.csv", index=False)
