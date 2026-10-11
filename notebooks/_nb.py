"""Tiny helpers for building notebooks programmatically."""

import nbformat as nbf


def notebook(cells):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    out = []
    for kind, text in cells:
        out.append(nbf.v4.new_markdown_cell(text.strip()) if kind == "md" else nbf.v4.new_code_cell(text.strip()))
    nb.cells = out
    return nb


def md(t):
    return ("md", t)


def code(t):
    return ("code", t)
