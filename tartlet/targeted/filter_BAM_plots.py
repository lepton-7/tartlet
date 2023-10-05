import click
import pickle
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.parsefuncs import plot_gen, is_interesting, bin_counts, gen_kernel


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
def main(pick_root, out_dir, bin_size):
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

    pass_rate_arr = comm.gather(pass_rate_local, root=0)

    if rank == 0:
        pass_rates = defaultdict(list)
        for instance_dict in pass_rate_arr:
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
