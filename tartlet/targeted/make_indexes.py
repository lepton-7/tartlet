import click

from glob import glob
from pathlib import Path
from subprocess import run, PIPE
from tartlet.utils.utils import print


@click.command()
@click.option(
    "-i",
    "--ref-dir",
    required=True,
    help="Directory with fasta files to build indexes. Fasta files must have the .fna extension.",
)
@click.option(
    "-p",
    "--threads",
    default=4,
    show_default=True,
    help="Number of hisat2-build threads to launch.",
)
def main(ref_dir, threads):
    arr = glob(f"{ref_dir}/*.fna")

    print(f"Found these reference files:")
    for i in arr:
        print(i)

    print("-------------------------------------------------")
    print("-------------------------------------------------\n")

    for fasta in arr:
        fasta = Path(fasta)
        fasta_name = fasta.stem
        index_dir = fasta.parent.joinpath(f"{fasta.stem}_index")

        Path(index_dir).mkdir(exist_ok=True, parents=True)

        print(f"Generating index for {fasta}")

        call = run(
            [
                "hisat2-build",
                "-p",
                f"{threads}",
                f"{fasta}",
                f"{index_dir}/{fasta_name}_index",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        if call.returncode:
            print(f"Failed:\n{call.stderr.decode('utf-8')}")
        else:
            print("Success")
            print(call.stdout.decode("utf-8"))
