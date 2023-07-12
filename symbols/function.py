from typing import Tuple


class Expression:

    def __init__(self, args: Tuple) -> None:
        self._args = args

    def __str__(self):
        return f'{self._args[0]} {self._args[1]} {self._args[2]}'