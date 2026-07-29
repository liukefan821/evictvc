# Revision working copy

Derived from `../submitted-cic2026_3/evictvc_iacrcc_submission.tex`.
Still `\documentclass[version=submission]{iacrcc}`; the submitted PDF is
reproducible from the frozen directory, and this one differs from it in
the bibliography only.

## What changed

The hand-written `thebibliography` block (11 `\bibitem` entries) was
replaced by `\bibliography{abbrev3,crypto,biblio}`.

- `abbrev3.bib`, `crypto.bib` — CryptoBib, copied from
  github.com/cryptobib/export at commit 84bfd10. Both are needed:
  `crypto.bib` uses `@string` macros defined in `abbrev3.bib`.
- `biblio.bib` — the seven ML references CryptoBib does not carry.

Five citations now use CryptoBib keys: `AC:KatZavGol10`,
`EPRINT:FeiKho23` (the ePrint 2023/033 version, which carries a bug fix
over the 2020 note), `PKC:CatFio13`, `SCN:TABDFK20`, and `EC:BonBoy04a`
(newly cited at Assumption 6.1, which previously stated q-SDH without
attribution). The seven `biblio.bib` entries keep their original keys, so
the corresponding `\cite` calls are unchanged.

`iacrcc.cls` forces `\bibliographystyle{alphaurl}` and redefines
`\bibliographystyle` to raise a ClassError, so labels are alphabetic
([KZG10], [BB04]) rather than numeric, and DOIs are printed as links.
Result: 14 pages, 12 entries.

## Open items

- TopLoc: page range within PMLR 267 not established; the `pages` field
  is absent.
- VeriLLM: `author = {Ke Wang and others}` produces the label [W+25].
- `\begin{abstract}` (line 81) and `\begin{textabstract}` (line 104) do
  not match. `textabstract` is the field IACR extracts for metadata; a
  third variant was pasted into the CiC submission form, with PDF line
  numbers embedded in it.
- `.tex` line 803 carries a `% TODO` comment about the Section 9 prose.

## Build

    ./build.sh
