import click
import importlib_resources
import pathlib


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
