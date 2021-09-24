"""Test the :mod:`sphinx` documentation generation."""

import pathlib

import pytest
from sphinx.application import Sphinx
from sphinx.util.docutils import docutils_namespace


@pytest.mark.parametrize("buildername", ["html", "latex", "epub"])
def test_sphinx_build(buildername: str, tmp_path: pathlib.Path) -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    doctreedir = tmp_path / "doctrees"
    with docutils_namespace():
        app = Sphinx("docs", "docs", tmp_path, doctreedir,
                    buildername=buildername, warningiserror=True)
        app.build(force_all=True)
