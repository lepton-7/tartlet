import click

from glob import glob
from pathlib import Path
from subprocess import run, PIPE
from tart.utils.utils import print


@click.command(context_settings={"ignore_unknown_options": True})
@click.option(
    "-i",
    "--ref-dir",
    required=True,
    help="Directory with fasta files and built indexes. Fasta files must have the .fna extension.",
)
@click.option("-1", "--m1", required=True, help="Path to the first-in-mate reads.")
@click.option("-2", "--m2", required=True, help="Path to the second-in-mate reads.")
@click.option(
    "-o",
    "--out-dir",
    required=True,
    help="Directory to place SAM output. SAM files will be further organised into subdirectories.",
)
@click.option(
    "--readset-name",
    default=None,
    show_default=True,
    help="Read pair that is added to SAM filename. If left undefined, will be the first-in-mate filename",
)
@click.argument(
    "hisat2",
    nargs=-1,
)
def main(ref_dir, m1, m2, out_dir, readpair_name, hisat2):
    """Pass options to the HISAT2 invocation."""
    refs = glob(f"{ref_dir}/*.fna")

    if readpair_name is None:
        # Get the name of the first mate file without extension
        readpair_name = Path(m1).name.split(".")[0]

    for fasta in refs:
        idx_dir = f"{fasta[:-4]}_index"
        rswtch_class = Path(fasta).name[:-4]

        samout_dir = f"{out_dir}/{rswtch_class}"
        Path(samout_dir).mkdir(parents=True, exist_ok=True)

        samout_path = f"{samout_dir}/{readpair_name}.sam"

        call = run(
            [
                "hisat2",
                "-x",
                f"{idx_dir}/{rswtch_class}_index",
                "-1",
                m1,
                "-2",
                m2,
                "-S",
                samout_path,
                *hisat2,  # options fed through from function call
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        if call.returncode:
            print(
                f"Failed:\n{call.stderr.decode('utf-8')}",
            )

        else:
            print(call.stdout.decode("utf-8"))

        print(
            "\n---------------------------------------------------------------------\n"
        )

    print("Done")
