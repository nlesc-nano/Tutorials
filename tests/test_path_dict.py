from __future__ import annotations

import os
import sys
import weakref
import pickle
from types import MappingProxyType
from pathlib import Path
from typing import Tuple, Any, TYPE_CHECKING, Generator, TypeVar, Callable, Union, Type

import pytest
from tutorials import PathDict

if TYPE_CHECKING:
    from _pytest.mark.structures import ParameterSet

T = TypeVar("T")
_PathLike = Union[str, bytes, os.PathLike]

DCT = PathDict({"a": 1, "b": 2, "c": 3})
DCT_REF = PathDict({"a": 1, "b": 2, "c": 3, "d": 4})
DCT_REF2 = PathDict({"a": None, "b": None, "c": None})
DCT_REF3 = PathDict({"a": 1, "b": 1, "c": 1})

LT_39 = pytest.mark.skipif(sys.version_info < (3, 9), reason="requires python 3.9")


def get_param(
    n: None | int = None,
    **kwargs: Any | Tuple[Any, ...],
) -> Generator[ParameterSet, None, None]:
    """Like normal :func:`pytest.mark.parametrize`, except the first value field is used as ``id``."""
    for name, args in kwargs.items():
        if not isinstance(args, tuple):
            args = (args,)
        if n is not None and len(args) == n + 1:
            *args, marks = args
            yield pytest.param(*args, marks=marks, id=name)
        else:
            yield pytest.param(*args, id=name)


def _test_setitem(dct: PathDict[T], key: _PathLike, value: T) -> T:
    __tracebackhide__ = True

    dct[key] = value
    return dct[key]


def _test_delitem(dct: PathDict[T], key: _PathLike) -> bool:
    __tracebackhide__ = True

    del dct[key]
    try:
        dct[key]
    except KeyError:
        return True
    else:
        return False


def _test_update(dct: PathDict[T], value: Any) -> PathDict[T]:
    __tracebackhide__ = True

    dct.update(value)
    return dct


def _test_ior(dct: PathDict[T], value: Any) -> PathDict[T]:
    __tracebackhide__ = True

    dct |= value
    return dct


class TestPathDict:
    """Tests for the :class:`tutorials.PathDict` class."""

    @pytest.mark.parametrize("func,ref", get_param(
        n=2,
        as_pathlib=(lambda dct, func: dct.as_pathlib(), {Path(k): v for k, v in DCT.items()}),
        setdefault1=(lambda dct, func: dct.setdefault(func("a")), 1),
        setdefault2=(lambda dct, func: dct.setdefault(func("d")), None),
        setdefault3=(lambda dct, func: dct.setdefault(func("d"), 2), 2),
        get1=(lambda dct, func: dct.get(func("a")), 1),
        get2=(lambda dct, func: dct.get(func("d")), None),
        get3=(lambda dct, func: dct.get(func("d"), 2), 2),
        pop1=(lambda dct, func: dct.pop(func("a")), 1),
        pop2=(lambda dct, func: dct.pop(func("d"), 2), 2),
        update1=(lambda dct, func: _test_update(dct, {func("d"): 4}), DCT_REF),
        update2=(lambda dct, func: _test_update(dct, [(func("d"), 4)]), DCT_REF),
        copy=(lambda dct, func: dct.copy(), DCT),
        pickle=(lambda dct, func: pickle.loads(pickle.dumps(dct)), DCT),
        weakref=(lambda dct, func: weakref.ref(dct)(), DCT),
        fromkeys1=(lambda dct, func: dct.fromkeys(DCT.keys()), DCT_REF2),
        fromkeys2=(lambda dct, func: dct.fromkeys(DCT.keys(), 1), DCT_REF3),
        __getitem__=(lambda dct, func: dct[func("a")], 1),
        __setitem__=(lambda dct, func: _test_setitem(dct, func("a"), 2), 2),
        __delitem__=(lambda dct, func: _test_delitem(dct, func("a")), True),
        __contains__1=(lambda dct, func: func("a") in dct, True),
        __contains__2=(lambda dct, func: func("d") in dct, False),
        __contains__3=(lambda dct, func: 1 in dct, False),
        __or__1=(lambda dct, func: dct | {func("d"): 4}, DCT_REF, LT_39),
        __or__2=(lambda dct, func: dct | [(func("d"), 4)], DCT_REF, LT_39),
        __ior__1=(lambda dct, func: _test_ior(dct, {func("d"): 4}), DCT_REF, LT_39),
        __ior__2=(lambda dct, func: _test_ior(dct, [(func("d"), 4)]), DCT_REF, LT_39),
        __init__1=(lambda dct, func: type(dct)({"a": 1, "b": 2, "c": 3}), DCT),
        __init__2=(lambda dct, func: type(dct)(a=1, b=2, c=3), DCT),
        __init__3=(lambda dct, func: type(dct)([("a", 1), ("b", 2), ("c", 3)]), DCT),
        __repr__1=(lambda dct, func: repr(dct).__class__, str),
        __repr__2=(lambda dct, func: repr(PathDict({str(i): i for i in range(99)})).__class__, str),
        __eq__1=(lambda dct, func: {func(k): v for k, v in dct.items()} == DCT, True),
        __eq__2=(lambda dct, func: MappingProxyType({func(k): v for k, v in dct.items()}) == DCT, True),
        __ne__1=(lambda dct, func: dct != 1, True),
        __ne__2=(lambda dct, func: dct != {1: 1}, True),
    ))
    @pytest.mark.parametrize("key_converter", get_param(
        str=lambda n: n,
        bytes=lambda n: n.encode(),
        Path=lambda n: Path(n),
    ))
    def test_pass(
        self,
        func: Callable[[PathDict[Any], Callable[[str], _PathLike]], T],
        ref: T,
        key_converter: Callable[[str], _PathLike],
    ) -> None:
        """Test whether the supplied operation passes as expected."""
        dct = DCT.copy()
        output = func(dct, key_converter)
        assert output == ref

    @pytest.mark.parametrize("func,exc", get_param(
        n=2,
        __init__1=(lambda: PathDict({1: 1}), TypeError),
        __init__2=(lambda: PathDict(1), TypeError),
        setdefault=(lambda: DCT.setdefault(1), TypeError),
        get=(lambda: DCT.get(1), TypeError),
        pop=(lambda: DCT.pop(1), TypeError),
        update1=(lambda: DCT.update({1: 1}), TypeError),
        update2=(lambda: DCT.update(1), TypeError),
        fromkeys=(lambda: PathDict.fromkeys([1]), TypeError),
        __getitem__=(lambda: DCT[1], TypeError),
        __setitem__=(lambda: _test_setitem(DCT, 1, None), TypeError),
        __delitem__=(lambda: _test_delitem(DCT, 1), TypeError),
        __or__1=(lambda: DCT | {1: 1}, TypeError, LT_39),
        __or__2=(lambda: DCT | 1, TypeError, LT_39),
        __ior__1=(lambda: _test_ior(DCT, {1: 1}), TypeError, LT_39),
        __ior__2=(lambda: _test_ior(DCT, 1), TypeError, LT_39),
    ))
    def test_raises(self, func: Callable[[], Any], exc: Type[Exception]) -> None:
        """Test whether the supplied operation raises the expected error."""
        with pytest.raises(exc):
            func()
