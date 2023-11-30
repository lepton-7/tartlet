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

from tart.utils import (
    infernal,
    prodigal,
)


@click.group()
def targeted():
    pass


@click.group()
def utils():
    pass


targeted.add_command(align_reads.main, name="align")
targeted.add_command(filter_BAM_plots.exec_main, name="filter")
targeted.add_command(make_BAMs.main, name="convert-sam")
targeted.add_command(make_indexes.main, name="index")
targeted.add_command(make_reference_sequences.main, name="reference-gen")
targeted.add_command(parse_BAMs.main, name="parse-bam")
targeted.add_command(switch_loc_in_ref.main, name="bounds")

utils.add_command(infernal.default_scan_for_riboswitches, name="find-riboswitches")
utils.add_command(prodigal.default_prodigal, name="find-orfs")
utils.add_command(prodigal.record_orf_locations, name="record-downstream")
