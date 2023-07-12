from typing import Any
from .function import Expression


class Args:

    def __init__(self) -> None:
        pass


class Symbol:

    def __init__(self, arg: str) -> None:
        self._symbol = arg

    def __pow__(self, arg) -> str:
        return Expression((self._symbol, '**', arg))

    def __call__(self) -> Any:
        return self._symbol
    
    def __str__(self) -> str:
        return self._symbol


