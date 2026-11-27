# EduMathPhys

Japanese self-study site for university mathematics and physics (undergraduate to graduate), built with Astro + Starlight.

## Commands

- `npm run dev` — dev server
- `npm run check` — fast consistency check (registry, labels, refs, citations, terms, KaTeX)
- `npm run build` — full build (fails on any content error, KaTeX error or broken link)
- `npm run typecheck` — `astro check`

## Layout

- `curriculum/subjects/<area>/<subject>.yaml` — the registry: subjects → chapters → sections, prerequisites, `status: stub|draft|reviewed`. Source of truth for navigation, stub pages, roadmap and coverage.
- `curriculum/references.yaml` — bibliography (well-known textbooks only).
- `src/content/docs/<area>/<subject>/<chapter-id>.md` — chapter text. A file must exist iff status is not `stub`.
- `src/plugins/remark-edumath.ts` — directives (theorem envs, refs, equations, exercises, quizzes, citations, terms, TikZ, interactive).
- `src/lib/` — curriculum loader, document index (labels/terms/notation), sidebar, graph layout, KaTeX macros.
- `src/pages/` — generated pages (home, subject index, stub chapters, map, glossary, notation, coverage, progress).
- `docs/STYLE.md` — writing rules and markup reference. Read it before writing content.

## Writing content

Follow `docs/STYLE.md`. Writing priority (continue from the first subject with stub chapters):

1. linear-algebra, analysis, analytical-mechanics, quantum-mechanics
2. topology, electromagnetism, group-theory, thermodynamics, complex-analysis, statistical-mechanics,
   measure-integration, special-relativity, ode, advanced-quantum-mechanics, rings-modules, fields-galois,
   functional-analysis, mathematical-physics, manifolds, general-relativity, probability, quantum-field-theory, solid-state
3. all remaining subjects

Workflow: chapters are drafted in parallel by Codex subagents (one subject per agent, a few chapters per batch, `status: draft`), then reviewed by a Claude reviewer agent that verifies proofs and computations (`uv run --with sympy python ...`) and fixes errors. After review, run `npm run check && npm run build`, set `status: reviewed`, and commit only the reviewed files.

## Rules

- Comments, commit messages and docs for developers in English; site content in Japanese.
- Generated artifacts (`dist/`, `.astro/`, `.cache/`) are gitignored.
- Never hard-code theorem/equation numbers in text; always use `:ref` / `:eqref`.
- Never guess bibliographic details.
