from __future__ import annotations

import os
import sys
import types
import pprint
import reprlib
import textwrap
import pathlib
from typing import Dict, TypeVar, ClassVar, Callable, Type, Any, Mapping
from collections.abc import Iterable

if sys.version_info >= (3, 8, 0):
    import functools
    pformat: Callable[[object], str] = functools.partial(pprint.pformat, sort_dicts=False)
else:
    pformat = pprint.pformat

_VT = TypeVar("_VT")
_TT = TypeVar("_TT", bound=Type[Any])

__all__ = ["PathDict"]


class _NoValue:
    __slots__ = ("__weakref__",)
    _cache: ClassVar[_NoValue] = NotImplemented

    def __new__(cls) -> _NoValue:
        if cls._cache is NotImplemented:
            cls._cache = object.__new__(cls)
        return cls._cache

    def __repr__(self) -> str:
        return "<no value>"


_NO_VALUE = _NoValue()


def _copy_docstring(cls1: Type[Any]) -> Callable[[_TT], Callable[[], _TT]]:
    """A class-decorator for copying all embedded docstrings from **cls1** into **cls2**."""
    def set_docstring(cls2: _TT) -> Callable[[], _TT]:
        attr_intersection = vars(cls1).keys() & vars(cls2).keys()
        for name in attr_intersection:
            attr1: Callable[..., Any] = getattr(cls1, name)
            attr2: Callable[..., Any] = getattr(cls2, name)
            doc1: None | str = getattr(attr1, "__doc__", None)
            doc2: None | str = getattr(attr2, "__doc__", None)
            if doc1 is not None and doc2 is None:
                if isinstance(attr2, types.MethodType):
                    attr2.__func__.__doc__ = attr1.__doc__
                else:
                    attr2.__doc__ = attr1.__doc__
        return cls2
    return set_docstring


def _is_mapping_like(obj):
    """Check if **obj** has the :meth:`~dict.keys` and :meth:`~dict.__getitem__` methods."""
    return (
        callable(getattr(obj, "keys", None)) and
        callable(getattr(obj, "__getitem__", None))
    )


@_copy_docstring(dict)
class PathDict(Dict[str, _VT]):
    """A :class:`dict` subclass that takes :term:`path-like object <path-like objects>` as keys and automatically converts them into strings."""

    __slots__ = ("__weakref__",)

    def __init__(self, iterable=None, **kwargs) -> None:
        if iterable is None:
            return super().__init__(**kwargs)

        if _is_mapping_like(iterable):
            iterator = ((os.fsdecode(k), iterable[k]) for k in iterable.keys())
        elif isinstance(iterable, Iterable):
            iterator = ((os.fsdecode(k), v) for k, v in iterable)
        else:
            raise TypeError(f"{iterable.__class__.__name__!r} object is not iterable")
        return super().__init__(iterator, **kwargs)

    def as_pathlib(self) -> Dict[pathlib.Path, _VT]:
        """Return a new dictionary wherein all keys are converted into :class:`pathlib.Path` instances."""
        return {pathlib.Path(k): v for k, v in self.items()}

    def setdefault(self, key, default=None):
        key = os.fsdecode(key)
        return super().setdefault(key, default)

    def get(self, k, default=_NO_VALUE):
        k = os.fsdecode(k)
        return super().get(k) if default is _NO_VALUE else super().get(k, default)

    def pop(self, k, default=_NO_VALUE):
        k = os.fsdecode(k)
        return super().pop(k) if default is _NO_VALUE else super().pop(k, default)

    def update(self, iterable=None, **kwargs) -> None:
        if iterable is None:
            return super().update(**kwargs)

        if _is_mapping_like(iterable):
            iterator = ((os.fsdecode(k), iterable[k]) for k in iterable.keys())
        elif isinstance(iterable, Iterable):
            iterator = ((os.fsdecode(k), v) for k, v in iterable)
        else:
            raise TypeError(f"{iterable.__class__.__name__!r} object is not iterable")
        return super().update(iterator, **kwargs)

    def copy(self):
        # Bypass the relativelly slow `PathDict` constructor
        cls = type(self)
        ret = dict.__new__(cls)
        dict.__init__(ret, self)
        return ret

    @classmethod
    def fromkeys(cls, iterable, value=None):
        iterator = (os.fsdecode(i) for i in iterable)
        return super().fromkeys(iterator, value)

    def __getitem__(self, k):
        k = os.fsdecode(k)
        return super().__getitem__(k)

    def __setitem__(self, k, v):
        k = os.fsdecode(k)
        return super().__setitem__(k, v)

    def __delitem__(self, k):
        k = os.fsdecode(k)
        return super().__delitem__(k)

    def __contains__(self, k):
        try:
            k = os.fsdecode(k)
        except TypeError:
            pass
        return super().__contains__(k)

    @reprlib.recursive_repr(fillvalue="...")
    def __repr__(self):
        ret = pformat(dict(self))
        cls_name = self.__class__.__name__
        if "\n" in ret:
            indent = " " * (len(cls_name) + 1)
            return f"{cls_name}({textwrap.indent(ret, indent)[len(indent):]})"
        else:
            return f"{cls_name}({ret})"

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, (dict, types.MappingProxyType)):
            return NotImplemented

        cls = type(self)
        try:
            value_parsed: Mapping[Any, Any] = cls(value) if not isinstance(value, cls) else value
        except Exception:
            value_parsed = value
        return super().__eq__(value_parsed)

    if sys.version_info >= (3, 9):
        def __or__(self, value):
            ret = self.copy()
            ret.update(value)
            return ret
