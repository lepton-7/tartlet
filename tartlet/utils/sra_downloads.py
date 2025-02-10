import click

from pathlib import Path
from subprocess import run, PIPE
from tartlet.utils.utils import print
from tartlet.utils.mpi_context import BasicMPIContext


@click.command()
@click.option(
    "-i",
    "--acc-file",
    required=True,
    help="List of SRA accession IDs, one on each line.",
)
@click.option(
    "-t",
    "--temp-dir",
    required=True,
    help="Temporary directory to use when downloading un-compressed sequence data.",
)
@click.option(
    "-o",
    "--save-dir",
    required=True,
    help="Output directory for the fastq.gz files.",
)
def main(acc_file, temp_dir, save_dir):

    # Make sure dirs exist
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    Path(temp_dir).mkdir(parents=True, exist_ok=True)

    # Determine MPI context ---------------------------------------------------
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

    if not Path(acc_file).exists() or Path(acc_file).is_dir():
        if rank == 0:
            print(f"Input list path ({acc_file}) is invalid, exiting.")

        return

    # Read in the list
    with open(acc_file, "r") as f:
        acc_list = [x for x in f.read().splitlines()]

    if len(acc_list) == 0:
        if rank == 0:
            print(f"Input list {acc_file} is empty, exiting.")

        return

    mp_con.set_full_list(acc_list)
    worker_list: list[str] = mp_con.generate_worker_list()

    maxTries = 3
    for sra in worker_list:
        i = 0
        failed = True
        while i < maxTries and failed:
            c1 = run(
                [
                    "prefetch",
                    f"{sra}",
                    "-O",
                    f"{temp_dir}/SRAs/",
                ],
                stdout=PIPE,
                stderr=PIPE,
            )
            failed = False

            if c1.returncode:
                print(f"Retrying {sra} fetch ({i}); Error code ({c1.returncode})")
                print(c1.stderr)

                i += 1
                failed = True

        i = 0
        failed = True
        while i < maxTries and failed:
            c2 = run(
                [
                    "fasterq-dump",
                    f"{temp_dir}/SRAs/{sra}/{sra}.sra",
                    "--split-files",
                    "--temp",
                    f"{temp_dir}",
                    "--outdir",
                    f"{temp_dir}",
                ],
                stdout=PIPE,
                stderr=PIPE,
            )
            failed = False

            if c2.returncode:
                print(f"Retrying {sra} fastq dump ({i})")
                print(c2.stderr)

                i += 1
                failed = True

        i = 0
        failed = True
        while i < maxTries and failed:
            c3 = run(
                [
                    "gzip",
                    "--best",
                    f"{temp_dir}/{sra}_1.fastq",
                    f"{temp_dir}/{sra}_2.fastq",
                ],
                stdout=None,
                stderr=None,
            )
            failed = False

        try:  # Move the compressed archives back from temp
            movrun = run(
                [
                    "mv",
                    f"{temp_dir}/{sra}_1.fastq.gz",
                    f"{temp_dir}/{sra}_2.fastq.gz",
                    f"{save_dir}/",
                ],
                stdout=PIPE,
                stderr=PIPE,
            )
        except:
            print(f"mv failed for {sra}: \n{movrun.stderr}")

        try:
            delrun = run(
                [
                    "rm",
                    "-r",
                    f"{temp_dir}/SRAs/{sra}",
                ],
                stdout=PIPE,
                stderr=PIPE,
            )

        except:
            print(f"rm failed: \n{delrun.stderr}")

        print(f"Success: {sra}")
