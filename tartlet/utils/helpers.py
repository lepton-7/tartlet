import click


def print(obj):
    """Wrapper for click.echo to override default print.

    Args:
        obj (Any): Object to print
    """
    click.echo(obj)
