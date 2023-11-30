import click
import pathlib
import importlib_resources

from pandas import Series


def print(obj):
    """Wrapper for click.echo to override default print.

    Args:
        obj (Any): Object to print
    """
    click.echo(obj)


def get_datapath_obj() -> pathlib.Path:
    """Get a path for the data subdirectory independent of package installation knowledge.

    Returns:
        pathlib.Path: Data subdirectory
    """
    utils = importlib_resources.files("tart.utils")
    with importlib_resources.as_file(utils / "data") as path:
        return path


def rowid(row: Series) -> str:
    """Generate the rowid for a ledger row.

    This should be used across the entire package when determining the rowids.

    Args:
        row (Series): Row element for which to generate the ID

    Returns:
        str: rowid
    """
    classname = row["target_name"]
    qname = row["query_name"]
    frm = row["seq_from"]
    to = row["seq_to"]
    strand = row["strand"]

    # No spaces to ensure the entire string is recognised as the ID
    return f"{classname}#{qname}#{frm}#{to}#{strand}"
