import click

from tart.targeted import (
    align_reads,
    filter_BAM_plots,
    make_BAMs,
    make_indexes,
    make_reference_sequences,
    parse_BAMs,
    switch_loc_in_ref,
)


@click.group()
def cli():
    pass


cli.add_command(align_reads.main, name="align")
cli.add_command(filter_BAM_plots.main, name="filter")
cli.add_command(make_BAMs.main, name="convert-sam")
cli.add_command(make_indexes.main, name="index")
cli.add_command(make_reference_sequences.main, name="reference-gen")
cli.add_command(parse_BAMs.main, name="parse-bam")
cli.add_command(switch_loc_in_ref.main, name="bounds")
