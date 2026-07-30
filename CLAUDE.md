# CLAUDE.md



---

## 0. Project snapshot

<!-- Fill in / keep current -->
- **Domain:** population genomics / computational genomics (WGS pipelines, tree-sequence inference, ABC/Bayesian ML methods)
- **Languages:** R (packages, stats), Python (pipelines, ML), Snakemake (workflow orchestration), bash (glue)
- **Data:** VCFs, BAMs, tree sequences (`.trees`), large reference genomes — treat all of these as too big to `cat`/load blindly
- **Entry points:** `Snakefile`, `workflow/rules/*.smk`, `R/`, `scripts/`
- **Env managers:** `uv` (Python), `renv` or `conda` (R) — check which one is active before installing anything

---

## 1. How Claude should work here (general best practices)

These are drawn from Anthropic's own published Claude Code guidance — not project-specific quirks, so don't remove them casually.

- **Explore → Plan → Code → Commit.** For anything non-trivial, read relevant files and state a plan *before* editing. Use plan mode for multi-file or multi-step changes. Don't jump straight to code on ambiguous asks — ask one clarifying question if genuinely blocked, otherwise state the assumption and proceed.
- **Small, verifiable steps.** Prefer several small diffs with checks in between over one large diff. Run tests/linters after each meaningful change, not just at the end.
- **Be honest about uncertainty.** If a fix works "by construction" but you haven't actually run it, say so explicitly. Never claim a test passed, a script ran, or a number was computed without having actually executed it in this session.
- **Don't silently paper over failures.** If a command errors, a package is missing, or an assumption breaks, stop and report it — don't quietly work around it in a way that changes the science (e.g. dropping failing samples, changing a filter threshold) without flagging it.
- **Use subagents/parallel exploration for large search tasks** (e.g. "find every place ploidy is assumed to be 2") rather than trying to hold it all in one context window.
- **Keep context lean.** Use `/clear` between unrelated tasks. Don't paste entire large files into chat when a targeted `grep`/`view` of the relevant section will do.
- **Git hygiene:** small commits, descriptive messages, never `git push --force` or rewrite shared history without explicit confirmation. Never commit data files, credentials, or `.Renviron`/`.env` secrets — check `.gitignore` covers them.
- **Ask before anything destructive or irreversible:** deleting files, overwriting pipeline outputs, dropping database tables, force-pushing, running on the full dataset instead of a test subset.

---

## 2. Sanity guidelines for data science / stats work

- **Reproducibility first.** Every stochastic step (bootstrap, MCMC, train/test split, ABC simulation, neural net init) must have an explicit, logged seed. No "it worked when I ran it" results without a fixed seed to reproduce them.
- **Never fabricate numbers, citations, or results.** If a benchmark, p-value, or effect size isn't actually computed from real output in this session, don't present it as if it were. If asked for a citation, verify it's real (don't hallucinate DOIs/author lists) — flag if unsure rather than guessing.
- **Sanity-check before interpreting.** Before trusting a result (e.g. a negative Ne estimate, a suspiciously perfect R², an empty VCF after filtering), check for the boring explanation first: coordinate off-by-one, wrong ploidy, silent NA propagation, unit mismatch (bp vs Mb, generations vs years).
- **State assumptions and caveats explicitly** in any methods text, especially around exchangeability assumptions (conformal prediction), model misspecification, or benchmark fairness — don't let a clean-looking result imply more certainty than the method supports.
- **Distinguish "ran without error" from "correct."** A pipeline finishing is not validation. Prefer a concrete sanity check (known-truth simulation, sanity plot, summary stat sense-check) over just "no exceptions thrown."
- **No p-hacking / silent re-running until significant.** If a test is repeated with different filters/subsets to change the outcome, that must be visible and justified, not hidden.

---

## 3. Genomics-specific gotchas

- **Coordinate systems:** VCF/GFF are 1-based; BED and many Python/tskit APIs are 0-based. State which convention is in use at every boundary (parsing, writing output, plotting).
- **Reference genome / assembly version** must be pinned and logged per run — never assume it matches a previous run.
- **Ancestral state / polarization:** flag when a site's ancestral allele is inferred vs. known; don't treat inferred ancestral states as ground truth downstream (tsinfer/tsdate especially).
- **Missing data ≠ absence of variation.** Don't silently treat missing genotypes as reference/ancestral.
- **Tree sequence provenance:** when producing/modifying a `.trees` file, preserve or update the provenance table so it's traceable which script/parameters produced it.
- **Per-scaffold / per-chromosome pipelines:** confirm outputs are being concatenated/compared in a consistent scaffold order — a silent misalignment here corrupts everything downstream without erroring.

---

## 4. Commands

<!-- Fill in and keep current — Claude will run these directly -->
```bash
# Environment
uv sync                          # Python deps
Rscript -e 'renv::restore()'     # R deps (if using renv)

# Pipeline
snakemake -n                     # dry run first, always
snakemake --cores 4

# Tests / checks
Rscript -e 'devtools::test()'    # R package tests (e.g. abcneuralnet)
pytest

# Lint / format
Rscript -e 'styler::style_pkg()'
ruff check .
```

---

## 5. Code style

- R: tidyverse-ish style, `styler`-compatible, roxygen2 docs for exported functions.
- Python: type hints where practical, `ruff`-clean.
- Prefer explicit over clever. A slower, readable filtering step beats a one-liner nobody can audit six months from now in a manuscript revision.
- Any function implementing a published method (loss function, statistic, distance measure) should cite the source (paper/equation number) in a comment or docstring.

---

## 6. When in doubt

Ask. A wrong assumption compounds silently across a Snakemake DAG or an MCMC run in a way that's expensive to catch later — a 10-second clarifying question is cheaper than a corrupted result discovered at manuscript revision stage.