from __future__ import annotations

import os
import sys
import types
import pprint
import reprlib
import textwrap
import pathlib
from typing import Dict, TypeVar, ClassVar, Callable, Type, Any, Mapping
from inspect import signature, Parameter, Signature
from collections.abc import Iterable

if sys.version_info >= (3, 8, 0):
    import functools
    pformat: Callable[[object], str] = functools.partial(pprint.pformat, sort_dicts=False)
else:
    pformat = pprint.pformat

_VT = TypeVar("_VT")
_TT = TypeVar("_TT", bound=Type[Any])

__all__ = ["PathMapping"]


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


def _modify_init_signature(func):
    """Change ``__iterable`` into a positional-only parameter."""
    sgn = signature(func)
    prm = list(sgn.parameters.values())
    prm[0] = prm[0].replace(kind=Parameter.POSITIONAL_ONLY)
    prm[1] = prm[1].replace(name="iterable", kind=Parameter.POSITIONAL_ONLY)
    func.__signature__ = sgn.replace(parameters=prm)
    return func


def _init(self, iterable, kwargs):
    self._hash = None
    if iterable is None:
        self._dict = dict(**kwargs)
        return None

    if _is_mapping_like(iterable):
        iterator = ((os.fsdecode(k), iterable[k]) for k in iterable.keys())
    elif isinstance(iterable, Iterable):
        iterator = ((os.fsdecode(k), v) for k, v in iterable)
    else:
        iterator = iterable
    self._dict = dict(iterator, **kwargs)


@_copy_docstring(dict)
class PathMapping(Mapping[str, _VT]):
    """A :class:`~collections.abc.Mapping` that takes :term:`path-like object <path-like objects>` as keys and automatically converts them into strings."""

    __slots__ = ("__weakref__", "_dict", "_hash")
    __signature__ = Signature([
        Parameter("iterable", Parameter.POSITIONAL_ONLY, default=None),
        Parameter("kwargs", Parameter.VAR_KEYWORD),
    ])

    if sys.version_info >= (3, 8, 0):
        def __init__(self, iterable=None, /, **kwargs):
            return self._init(self, iterable, kwargs)
    else:
        @_modify_init_signature
        def __init__(self, __iterable=None, **kwargs):
            return self._init(self, __iterable, kwargs)

    @classmethod
    def _reconstruct(cls, dct, hsh=None):
        """Construct a new **cls** instance while bypassing :meth:`__init__`."""
        ret = object.__new__(cls)
        ret._dict = dct
        ret._hash = hsh
        return ret

    def to_pathlib_dict(self):
        """Return a new dictionary wherein all keys are converted into :class:`pathlib.Path` instances."""
        return {pathlib.Path(k): v for k, v in self.items()}

    def keys(self):
        return self._dict.keys()

    def items(self):
        return self._dict.items()

    def values(self):
        return self._dict.values()

    def get(self, key, default=None):
        key = os.fsdecode(key)
        return self._dict.get(key, default)

    def __hash__(self):
        if self._hash is None:
            self._hash: None | int = hash(frozenset(self.items()))
        return self._hash

    def __copy__(self):
        return self

    def __deepcopy__(self, memo=None):
        return self

    def __reduce__(self):
        cls = type(self)
        return cls._reconstruct, (self._dict, self._hash)

    def __getitem__(self, key):
        key = os.fsdecode(key)
        return self._dict[key]

    def __contains__(self, key):
        try:
            key = os.fsdecode(key)
        except TypeError:
            pass
        return key in self._dict

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    @reprlib.recursive_repr(fillvalue="...")
    def __repr__(self):
        ret = pformat(self._dict)
        cls_name = self.__class__.__name__
        if "\n" in ret:
            indent = " " * (len(cls_name) + 1)
            ret_indent = textwrap.indent(ret, indent)[len(indent):]
            return f"{cls_name}({ret_indent})"
        else:
            return f"{cls_name}({ret})"

    def __eq__(self, value):
        return self._dict == value

    if sys.version_info >= (3, 9):
        def __or__(self, value):
            cls = type(self)
            other = cls(value)
            return cls._reconstruct(self._dict | other._dict)
