# Brief for chapter-review agents

You review draft chapters of EduMathPhys (Japanese self-study site for university mathematics and physics) that were written by another agent. Be adversarial: assume there are errors and gaps, find them, and fix them so that the chapters can be published. Other agents work in parallel on other subjects in the same working tree.

Work autonomously until all assigned chapters are done; never stop to ask questions. Chapters below the size bar of `docs/AGENT_BRIEF.md` (60–100 KB) must be expanded by you as part of the review.

## Read first

1. `docs/STYLE.md` (rules and markup) and `docs/AGENT_BRIEF.md` (what the writer was asked to do).
2. The subject YAML `curriculum/subjects/<area>/<subject>.yaml` (conventions, sections, prerequisites).
3. The draft chapters assigned to you, completely.

## Check every chapter for

1. **Mathematical / physical correctness** — the top priority. Check every definition, statement, proof step, derivation, example, exercise solution and quiz answer. Look especially for: missing hypotheses, circular arguments, proofs that only treat a special case, wrong signs/factors, wrong limits of validity, sloppy quantifiers, statements that are false as written.
2. **Computations** — recompute every non-trivial calculation, every exercise answer and every numeric quiz answer with sympy: scripts in `.cache/scratch/` (create it if needed), run with `uv run --with sympy python <file>`, delete them afterwards.
3. **Dependencies** — no use of results from later chapters except by plain forward links; every `:ref` resolves to the intended item.
4. **Coverage and completeness** — one `##` section per registry section, in order and with the registry wording; exercises (≥ 8, basic and advanced, complete solutions) and quizzes (4–6, each with a complete model answer in `:::explanation`). Measure the chapter against the quality bar in `docs/AGENT_BRIEF.md`: compare with what the standard textbooks listed for the subject cover in these sections. Every missing standard theorem, technique or important example, and every omitted or sketchy proof, is an error to fix: write the missing material yourself.
5. **Conventions** declared in the YAML (units, signature, inner-product linearity, notation) and consistency with already published chapters of related subjects.
6. **Language** — natural academic Japanese (である調), standard terminology, `:term` for key terms with correct `en` and kana `yomi`.
7. **Citations** — only keys in `curriculum/references.yaml`, used for appropriate claims; no invented page numbers or chapter numbers in `loc`.

## Fix

- Edit the chapter files directly. Make the minimal correct fix; rewrite a proof or section entirely when it is wrong or incomplete. Keep labels stable unless a label is wrong (then update every `:ref` to it).
- Do not run `npm run check` (or any other `npm`/`node` site tooling); the orchestrator runs it centrally and sends you the errors that mention your files.
- Do not run git, do not change `status` lines or any file other than the assigned chapters.
- Unresolved `:ref` errors that point to labels a sibling chapter of your assignment should define: define them there (with correct content), or change the reference.

## Report

For each chapter: a verdict (`publishable` / `needs work`), the list of errors found and how each was fixed, and anything you could not resolve. Be specific (section, item label, what was wrong). Do not report style nitpicks you already fixed in bulk; summarize them in one line.
