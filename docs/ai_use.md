# Use of AI assistance

© 2026 Syed Nafiz Reja Hasib.

This statement is for editors, reviewers and readers of this repository.

**Tool.** Claude (Anthropic), used through Claude Code, November 2026.

**What it did.**
- Drafted and tested the code in `src/`, `study/`, `verify/` and `tests/` to the author's specification in `docs/scope.md`.
- Ran the analyses, made the figures and exports, and drafted the README, the methods log (`docs/decisions.md`) and the report.
- Looked up the open sources listed in `references.md` and checked the values taken from them.

**What it did not do.**
- It is not an author.
- It did not set the research question, the hypothesis or its threshold; the author fixed these in `docs/scope.md` before any results were computed.
- Modelling choices it made under the author's standing instruction are listed in `docs/decisions.md` and marked "author to confirm".

**How the work was checked.**
- Every number in the README and the report comes from a file in `results/` or `verify/out/` written by a script in this repository; `python reproduce.py` rebuilds them.
- The verification checks in `verify/` (benchmarks against published and analytical solutions) ran before the study results, and the test suite runs on every push.
- The author reviewed the outputs and takes full responsibility for the content.
