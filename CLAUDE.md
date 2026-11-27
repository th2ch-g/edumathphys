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

Workflow (Claude token use is kept minimal; Codex does the heavy work):
1. Codex subagents draft chapters in parallel (one subject per agent, `status: draft`, brief in `docs/AGENT_BRIEF.md`); drafts under ~50 KB get a Codex expansion pass.
2. A separate Codex instance reviews and fixes them adversarially (`docs/REVIEW_BRIEF.md`).
3. Claude acts only as the final gate: `node scripts/stage-chapters.mjs <subject>/<chapter> ...` (check + build in a clean worktree with only those chapters), a light spot check, then commit and push.

## Rules

- Comments, commit messages and docs for developers in English; site content in Japanese.
- Generated artifacts (`dist/`, `.astro/`, `.cache/`) are gitignored.
- Never hard-code theorem/equation numbers in text; always use `:ref` / `:eqref`.
- Never guess bibliographic details.

## Resuming work (state at the last stop)

Unreviewed drafts live on the `drafts` branch (not published). To resume: `git checkout drafts -- src/content/docs curriculum/subjects` on top of main, then continue the workflow above.

- Written and expanded, waiting for Codex review: analysis ch4–13, analytical-mechanics ch4–12, quantum-mechanics (all 13).
- Written, not yet expanded or reviewed: complex-analysis, sets-logic, special-relativity, statistical-mechanics.
- Reviewed but held: electromagnetism/relativistic-em (waits for special-relativity to be published).
- Linear algebra ch1–2 predate the current quality bar (ch2 has only 6 exercises): needs an upgrade pass.
- References softened to plain links by `scripts/soften-refs.mjs` (restore as `:ref` once the target is published): electromagnetism/electrostatics (thm-anl-divergence), mechanics/kinematics, measure-integration/lebesgue-measure, measure-integration/differentiation. Find them with `grep -rn '](/math/\|](/physics/' src/content/docs` and compare with the labels in published chapters.
