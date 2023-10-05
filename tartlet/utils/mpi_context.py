import click
from typing import Any

try:
    from mpi4py import MPI

    no_mpi4py = False

except ImportError:
    click.echo("Mpi4py not found. Starting single-instance context")
    no_mpi4py = True


class PseudoMPICOMM:
    def __init__(self) -> None:
        pass

    def gather(self, sendobj: Any, root: int = 0):
        if root == 0:
            return [sendobj]
        else:
            raise ValueError("Root must be 0 for mimic gather")

    def bcast(self, obj: Any, root: int = 0):
        if root == 0:
            return obj
        else:
            raise ValueError("Root must be 0 for mimic bcast")


class BasicMPIContext:
    def __init__(self, full_list: list = None) -> None:
        self.full_list = full_list

        if no_mpi4py:
            self.is_SP = True
            self.comm = PseudoMPICOMM()
            self.size = 1
            self.rank = 0

        else:
            self.comm = MPI.COMM_WORLD
            self.size = self.comm.Get_size()
            self.rank = self.comm.Get_rank()

            self.is_SP = False if self.size > 1 else True

    def set_full_list(self, full_list: list) -> None:
        self.full_list = full_list

    def generate_worker_list(self) -> list:
        maxlen = len(self.full_list)

        # Determine the subset processed by one instance
        count = maxlen // self.size
        rem = maxlen % self.size

        if self.rank < rem:
            start = self.rank * (count + 1)
            stop = start + (count + 1)
        else:
            start = self.rank * count + rem
            stop = start + count

        self.worker_list = self.full_list[start:stop]

        # To indicate whether this worker instance is unecessary
        self.is_active = bool(len(self.worker_list))

        return self.worker_list
