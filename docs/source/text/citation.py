#!/usr/bin/env python3
"""Writes a BibTeX entry for libNeST from CITATION.cff.

Run by docs/source/Makefile as ``python text/citation.py CITATION.cff OUTPUT.bib``.
"""
import sys

import yaml


def bibtex(cff):
    authors = " and ".join("{}, {}".format(a["family-names"], a["given-names"])
                           for a in cff["authors"])
    fields = [("author", authors), ("title", cff["title"]), ("version", cff.get("version")),
              ("year", str(cff["date-released"])[:4] if "date-released" in cff else None),
              ("license", cff.get("license")), ("url", cff.get("repository-code"))]
    body = ",\n".join("  {:<8}= {{{}}}".format(k, v) for k, v in fields if v)
    return "@software{{libnest,\n{}\n}}\n".format(body)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: {} CITATION.cff OUTPUT.bib".format(sys.argv[0]))
    with open(sys.argv[1], encoding="utf-8") as f:
        cff = yaml.safe_load(f)
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        f.write(bibtex(cff))
