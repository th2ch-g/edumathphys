# Brief for chapter-writing agents

You are writing chapters for EduMathPhys, a Japanese self-study website for university mathematics and physics (Astro + Starlight). Several writers work in parallel in the same working tree, each on a different subject. Your assignment (subject and chapters) is given in your task message.

## Quality bar

The site must make standard university and graduate textbooks unnecessary: a reader who studies a subject here should not need to consult the books listed in `curriculum/references.yaml`. Concretely:

- Every definition, theorem and technique that the standard textbooks of the subject treat in the corresponding sections must appear, with complete proofs or derivations. Do not omit proofs as "standard" or "left to the reader"; the only acceptable omissions are results from a different subject (link to it) or results whose proof is genuinely far beyond the level (state it precisely, explain why the proof is omitted, and cite).
- Explain motivation and intuition before formal statements, add remarks on typical pitfalls and on how results are used later, and give several worked examples per section.
- Length follows content: typically 800–1500 lines per chapter. Split long proofs into lemmas.
- At least 8 exercises per chapter (basic and advanced, including proof problems and computations), all with complete solutions.

Work autonomously until the assignment is complete. Never stop to ask questions: when instructions conflict, the registry YAML wins (e.g. chapter titles exactly as written there, including spaces), then this brief, then STYLE.md; otherwise make the most reasonable decision and mention it in the final report.

## Read first

1. `docs/STYLE.md` — writing rules and the complete markup reference. Follow it exactly.
2. `src/content/docs/math/linear-algebra/linear-equations.md` and `vector-spaces.md` — finished reference chapters. Match their rigor, tone and markup; go beyond them in depth per the quality bar above.
3. `curriculum/subjects/<area>/<subject>.yaml` — your subject: conventions, chapter list, the section list of each chapter, prerequisites.
4. `curriculum/references.yaml` — the only allowed citation keys. Never add keys.
5. Existing chapters of prerequisite subjects (`src/content/docs/**`), to reuse their labels with `:ref` (list them with `grep -o '{#[a-z0-9-]*' <file>`).

## For each chapter

- File: `src/content/docs/<area>/<subject>/<chapter-id>.md`. Frontmatter `title` must equal the registry title exactly (including spaces); add `description` and a `notation` list.
- Structure: `## この章の目標` (3–5 bullets) → one `##` per registry section, in order and with the registry wording (`###` for subsections) → `## 演習問題` (at least 8 exercises, both `level=basic` and `level=advanced`, each with a complete worked `:::solution`, optional `:::hint`) → `## 確認テスト` (4–6 quizzes, mixing multiple choice and numeric `answer=`/`tol=`; each `:::explanation` is a complete model answer).
- Language: Japanese, である調, standard terminology; names of people in Latin script. Introduce key terms with `:term[用語]{en="english" yomi="ひらがな"}` (yomi: hiragana/katakana only).
- Rigor (mathematics): definitions before use; every theorem/proposition/lemma/corollary gets a complete proof. Cite with `:cite[key]` only for results genuinely beyond the chapter, and say why. Worked examples with concrete numbers.
- Rigor (physics): state assumptions and approximations; derivations complete and reproducible; units and orders of magnitude for physical quantities; mathematical statements get proofs.
- Labels: `{#<kind>-<subject prefix>-<name>}`, lowercase kebab-case, unique site-wide (kinds: def thm prop lem cor exm rem ax post law prin exr eq). Label every result or equation that is referenced later. Reference only labels that already exist (earlier chapters, other subjects' written chapters). Never hard-code theorem or equation numbers. For material that is not written yet, use plain links such as `[内積空間](/math/linear-algebra/inner-product-spaces/)`.
- Respect the conventions declared in the subject YAML (units, metric signature, inner-product linearity, ℏ explicit, etc.).
- KaTeX only: `\begin{aligned}` inside `$$` (no align/equation environments), no `\newcommand`, one `\label` per display block. Use the macros in STYLE.md.
- Figures: optional ```` ```tikz ```` / ```` ```tikzcd ```` code fences (plain TikZ; libraries: arrows.meta, calc, positioning, decorations.pathmorphing, decorations.markings, shapes.geometric). Keep them simple. Do not create interactive figures.
- Verify every non-trivial computation, every exercise answer and every numeric quiz answer with sympy: put throwaway scripts in `.cache/scratch/` and run `uv run --with sympy python .cache/scratch/<file>.py`. Delete your scripts when finished.
- In your subject YAML add the line `    status: draft` directly after the `title:` line of each chapter you created (before any `prerequisites:` line of that chapter). Change nothing else.
- Run `npm run check` after each chapter. It validates the whole site: fix every error that mentions your files and ignore errors in other subjects' files.

## Do not

- run `npm run build` or any git command;
- edit files other than your chapter files and the `status` lines of your subject YAML;
- add bibliography keys or create files under `src/` other than chapters.

## Final report

List the files created, the labels defined, the citations used, and every place where you are unsure about correctness.
