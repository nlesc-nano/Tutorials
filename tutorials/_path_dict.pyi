import os
import sys
import pathlib
from typing import TypeVar, Dict, Iterable, Tuple, overload, Union, Any
from typing_extensions import Protocol

_T = TypeVar("_T")
_KT = TypeVar("_KT")
_VT = TypeVar("_VT")
_VT_co = TypeVar("_VT_co", covariant=True)

_PathLike = Union[str, bytes, os.PathLike[Any]]

class _SupportsKeysAndGetItem(Protocol[_KT, _VT_co]):
    def keys(self) -> Iterable[_KT]: ...
    def __getitem__(self, __k: _KT) -> _VT_co: ...

__all__: list[str]

class PathDict(Dict[str, _VT]):
    @overload
    def __init__(self: PathDict[_VT]) -> None: ...
    @overload
    def __init__(self: PathDict[_VT], **kwargs: _VT) -> None: ...
    @overload
    def __init__(self, __map: _SupportsKeysAndGetItem[_PathLike, _VT], **kwargs: _VT) -> None: ...
    @overload
    def __init__(self, __iterable: Iterable[Tuple[_PathLike, _VT]], **kwargs: _VT) -> None: ...

    def copy(self) -> PathDict[_VT]: ...
    def setdefault(self, __key: _PathLike, __default: _VT = ...) -> _VT: ...
    def as_pathlib(self) -> Dict[pathlib.Path, _VT]: ...

    @overload
    def get(self, __key: _PathLike) -> None | _VT_co: ...
    @overload
    def get(self, __key: _PathLike, __default: _VT_co | _T) -> _VT_co | _T: ...

    @overload
    def pop(self, __key: _PathLike) -> _VT: ...
    @overload
    def pop(self, __key: _PathLike, __default: _VT_co | _T = ...) -> _VT_co | _T: ...

    @overload  # type: ignore[override]
    def update(self, __m: _SupportsKeysAndGetItem[_PathLike, _VT], **kwargs: _VT) -> None: ...
    @overload
    def update(self, __m: Iterable[Tuple[_PathLike, _VT]], **kwargs: _VT) -> None: ...
    @overload
    def update(self, **kwargs: _VT) -> None: ...

    @classmethod  # type: ignore[override]
    @overload
    def fromkeys(cls, __iterable: Iterable[_PathLike], __value: None = ...) -> PathDict[None | Any]: ...
    @classmethod
    @overload
    def fromkeys(cls, __iterable: Iterable[_PathLike], __value: _T) -> PathDict[_T]: ...

    def __getitem__(self, __k: _PathLike) -> _VT: ...
    def __setitem__(self, __k: _PathLike, __v: _VT) -> None: ...
    def __delitem__(self, __k: _PathLike) -> None: ...

    if sys.version_info >= (3, 9):
        def __or__(self, __value: _SupportsKeysAndGetItem[_PathLike, _T] | Iterable[Tuple[_PathLike, _T]]) -> PathDict[_VT | _T]: ...
        def __ior__(self, __value: _SupportsKeysAndGetItem[_PathLike, _VT] | Iterable[Tuple[_PathLike, _VT]]) -> PathDict[_VT]: ...
