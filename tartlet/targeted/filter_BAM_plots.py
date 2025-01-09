import click
import pickle
import pandas as pd
import tarfile

from glob import glob
from pathlib import Path
from tartlet.utils.cluster import Cluster
from tartlet.utils.plotting import CoveragePlot
from tartlet.utils.read_parsing import AlignDat, SegregatedAlignDat
from tartlet.utils.mpi_context import BasicMPIContext
from tartlet.utils.activity_inference import Candidate
from tartlet.utils.filter_functions import DefaultThresholds as checker
from tartlet.utils.activity_inference import get_peaks


def _log_cand_charac(peaklog: dict, cand: Candidate):
    peaklog["coverage_delta_absolute"] = cand.abs_cov_delta
    peaklog["coverage_delta_pval"] = cand.coverage_drop_pvalue
    peaklog["peak_summit"] = cand.summit
    peaklog["peak_width"] = cand.width
    # peaklog["peak_l/r"] = (cand.left, cand.right)
    # peaklog["coverage_delta_relative"] = cand.rel_cov_delta

    try:
        peaklog["coverage_delta_stable_relative"] = cand.stable_rel_cov_delta
    except AttributeError:
        pass

    peaklog["coverage_delta_noiseset"] = str(cand.coverage_delta_noise)

    peaklog["from_riboswitch_end"] = cand.from_switch_end
    peaklog["from_riboswitch_end_relative"] = cand.from_switch_end_relative

    # peaklog["ks_stat"] = cand.symks_stat
    # peaklog["ks_p"] = cand.symks_pval
    peaklog["decision_note"] = cand.note


def _process_peak(peak: Candidate, peaklog: dict, peaklog_loc: list):

    dec = checker.check(peak)
    peaklog["decision"] = dec

    _log_cand_charac(peaklog, peak)

    peaklog_loc.append(peaklog)

    return True if dec == "pass" else False


def _load_alignment_data(pick_root: str, pick_file: str, rank: int):
    """Helper to minimise nesting.

    Args:
        pick_root (str): Archive root.
        pick_file (str): Pickle filename to obtain from the archive.
        rank (int): MPI worker rank.

    Returns:
        segmented: SegregatedAlignDat object.
    """
    with tarfile.open(pick_root, "r:gz") as picktar:
        try:
            with picktar.extractfile(pick_file) as f:  # type: ignore
                segmented: SegregatedAlignDat = pickle.load(f)

                return segmented

        except KeyError:
            print(f"{pick_file} not found in archive {pick_root} on worker rank {rank}")

            return None


def _process_candidate_list(
    candlist: list[Candidate],
    align_charac: dict,
    charac_local: list,
    ref: str,
    transcriptome: str,
):
    """DEPR

    Args:
        candlist (list[Candidate]): _description_
        align_charac (dict): _description_
        charac_local (list): _description_
        ref (str): _description_
        transcriptome (str): _description_

    Returns:
        _type_: _description_
    """
    if len(candlist) == 0:
        align_charac["decision_note"] = "No suitable candidates"
        align_charac["decision"] = "fail"
        charac_local.append(align_charac)
        return "fail"

    for cand in candlist:
        align_charac = {}
        align_charac["rowid"] = ref
        align_charac["transcriptome"] = str(transcriptome)
        decision = checker.check(cand)
        align_charac["decision"] = decision

        if decision == "pass":
            charac_local.append(align_charac)
            _log_cand_charac(align_charac, cand)
            return decision

    # Getting to this point means there were no passes
    align_charac = {}
    align_charac["rowid"] = ref
    align_charac["transcriptome"] = str(transcriptome)
    decision = checker.check(candlist[0])
    align_charac["decision"] = decision
    _log_cand_charac(align_charac, candlist[0])
    charac_local.append(align_charac)

    return "fail"


@click.command()
@click.option(
    "-i",
    "--pick-root",
    required=True,
    help="Gunzipped archive for pickled outputs. Same file as the parser output.",
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
    default=15,
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
    "--noplots",
    is_flag=True,
    help="(Dev use) Suppresses all plot generation.",
)
@click.option(
    "--conv",
    is_flag=True,
    help="(Dev use) Output plots with the ends convolution panel.",
)
@click.option(
    "--statplot",
    is_flag=True,
    help="(Dev use) Output plot with stats.",
)
def exec_main(
    pick_root,
    out_dir,
    bin_size,
    min_cov_depth,
    ext_prop,
    run_depr,
    noplots,
    conv,
    statplot,
):
    if run_depr:
        # depr_main(pick_root, out_dir, bin_size, min_cov_depth, ext_prop, conv, statplot)
        raise ValueError("No deprecated function to run")

    else:
        main(
            pick_root,
            out_dir,
            bin_size,
            min_cov_depth,
            ext_prop,
            noplots,
            conv,
            statplot,
        )


def main(
    pick_root: str, out_dir, bin_size, min_cov_depth, ext_prop, noplots, conv, statplot
):
    # Determine MPI context ---------------------------------------------------
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

    # List of all pickled alignment data
    if rank == 0:
        try:
            with tarfile.open(pick_root, "r:gz") as picktar:
                total_files = [x.name for x in picktar.getmembers() if x.isfile()]
        except FileNotFoundError:
            print(f"Pickled data root {pick_root} not found.")
            total_files = None

    else:
        total_files = None

    total_files = comm.bcast(total_files, root=0)

    # Pickle root was not found
    if total_files is None:
        raise SystemExit(0)

    mp_con.set_full_list(total_files)
    worker_list: list[str] = mp_con.generate_worker_list()

    if rank == 0:
        print(f"Started {mp_con.size} instances.")
        print(f"Each instance running upto {len(worker_list)} iterations.")

    # Keep track of peak log for each alignment
    peak_log_local: list[dict] = []

    # Iterate through alignment datafiles -------------------------------------

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # First pass; find suitable candidates if present
    for pick_file in worker_list:
        # Load in the alignment data object from the archive
        segmented = _load_alignment_data(pick_root, pick_file, rank)
        if segmented is None:
            continue

        alignDat = segmented.total

        switch_size = alignDat.switch_end - alignDat.switch_start

        # Save path formation and sample/ref IDs
        ref = Path(pick_file[:-2]).name
        targetname = Path(pick_file).parts[-3]
        transcriptome = Path(pick_file).parts[-2]

        # Check whether there are enough reads to proceed.
        # This needs to be done to avoid clutter in the results
        # that failed the filter
        max_cov = max(alignDat.readcov[alignDat.switch_start : alignDat.switch_end])
        roi = (checker.relative_size_bounds, checker.relative_size_bounds)

        plotdec = "fail"

        # Iterate through each peak in the region of interest -----------------
        for peak in get_peaks(alignDat, ss_margins=ext_prop, reg_of_i=roi):
            # One row per peak
            peak_info = {}
            peak_info["rowid"] = ref
            peak_info["transcriptome"] = str(transcriptome)

            if max_cov < min_cov_depth:
                peak_info["decision_note"] = (
                    f"Minimum coverage not reached ({max_cov} < {min_cov_depth})"
                )
                peak_log_local.append(peak_info)
                break

            if _process_peak(peak, peak_info, peak_log_local):
                plotdec = "pass"

        if max_cov < min_cov_depth:
            continue

        passfaildir = plotdec

        # main plot path
        save_path = out_dir.joinpath(passfaildir, f"{ref}#{transcriptome}.png")
        save_path.parent.mkdir(exist_ok=True, parents=True)

        # stat plot path
        stat_path = out_dir.joinpath(passfaildir, f"{ref}#{transcriptome}_stats.png")

        # Calculate and set info for binned raw ends
        segmented.bin_rawends(bin_size=bin_size)

        # We still want to see a large region around the switch in the plots
        lplot = 1.0 if ext_prop[0] < 1.0 else ext_prop[0]
        rplot = 1.0 if ext_prop[1] < 1.0 else ext_prop[1]

        lbuff = int(switch_size * lplot)
        rbuff = int(switch_size * rplot)
        end_buffers = [lbuff, rbuff]

        pObj = CoveragePlot(segmented, end_buffers)
        if not noplots:
            pObj.default(save_path, with_conv=conv)

        if statplot and not noplots:
            pObj.distribution_plots(stat_path)

    if rank == 0:
        print("Waiting on all instances to start gather.")
    peak_log_arr = comm.gather(peak_log_local, root=0)

    if rank == 0:
        print("Completed gather.")
        if peak_log_arr is None:
            raise TypeError("Gather failed")
        log: list[dict] = []
        for instance_arr in peak_log_arr:
            instance_arr: list[dict]
            log.extend(instance_arr)

        # Make dataframe
        # df = pd.DataFrame({"target_name": classes, "pass_rate": rates})
        # df.to_csv(f"{out_dir}/pass_rates.csv", index=False)

        # Cluster peaks for later plotting
        peak_log, cluster_stats = Cluster(pd.DataFrame(log)).get()

        if len(peak_log) > 0:
            peak_log.to_csv(f"{out_dir}/peak_log.csv", index=False)
            cluster_stats.to_csv(f"{out_dir}/cluster_stats.csv", index=False)
        else:
            print(f"{pick_root} seems empty. No tables exported.")


def depr_main(
    pick_root: str, out_dir, bin_size, min_cov_depth, ext_prop, conv, statplot
):
    pass
