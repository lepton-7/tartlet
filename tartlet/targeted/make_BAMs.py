# Convert riboswitch alignment SAMs to sorted BAMs
import click

from sys import argv
from glob import glob
from subprocess import run, PIPE
from tart.utils.helpers import print
from tart.utils.mpi_context import BasicMPIContext


@click.command()
@click.option(
    "-i",
    "--sam-dir",
    required=True,
    help="Finds SAM files within its subdirectories and converts to sorted BAMS.",
)
def main(sam_dir):
    # List of all fastas in the input directory
    total_files = glob(f"{sam_dir}/**/*.sam")

    # MPI setup
    mp_con = BasicMPIContext(total_files)
    local_list = mp_con.generate_worker_list()

    # SAM -> sorted BAM commands
    for sam_file in local_list:
        file_base = sam_file[:-4]

        command_1_view = f"samtools view -S -b {sam_file}"
        command_2_sort = f"samtools sort {file_base}.bam -o {file_base}.sorted.bam"
        command_3_index = f"samtools index {file_base}.sorted.bam"

        with open(f"{file_base}.bam", "w") as outBAM:
            call_1 = run(command_1_view.split(" "), stdout=outBAM, stderr=PIPE)

        if call_1.returncode:
            print(call_1.stderr)

        call_2 = run(command_2_sort.split(" "), stdout=PIPE, stderr=PIPE)

        if call_2.returncode:
            print(call_2.stderr)

        call_3 = run(command_3_index.split(" "), stdout=PIPE, stderr=PIPE)

        if call_3.returncode:
            print(call_3.stderr)

        if call_1.returncode or call_2.returncode or call_3.returncode:
            print(f"Failed processing for: {sam_file}")
