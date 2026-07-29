#!/bin/sh
# Rebuild the revision PDF. Requires a TeX Live installation with
# pdflatex, latexmk and bibtex on PATH.
set -e
cd "$(dirname "$0")"
latexmk -g -pdf -pdflatex="pdflatex --no-shell-escape" evictvc_rev
