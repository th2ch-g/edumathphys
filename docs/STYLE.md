# Writing Style Guide

This guide defines how chapter content is written. All content is Japanese; code and comments are English.

## 1. Workflow per chapter

1. Read the chapter entry in `curriculum/subjects/<area>/<subject>.yaml` (title, sections, prerequisites).
2. Check existing labels in prerequisite chapters and reference them with `:ref[...]` instead of re-proving results.
3. Write `src/content/docs/<area>/<subject>/<chapter-id>.md`. The frontmatter `title` must equal the registry title.
4. Self-review with the checklist in section 7. Verify non-trivial computations with `uv run --with sympy python <script>` (scripts live in a scratch directory and are deleted afterwards).
5. Set `status: reviewed` (or `draft` if the checklist is not fully satisfied) in the registry.
6. Run `npm run check && npm run build`.

## 2. Chapter layout

```markdown
---
title: <registry title>
description: <one sentence summary>
notation:
  - symbol: '\Ker f'
    meaning: 線形写像 f の核
---

## この章の目標

- 箇条書きで 3〜5 項目

## <section 1 from the registry>
...

## 演習問題

## 確認テスト
```

- Sections (`##`) follow the registry `sections` list in order. Subsections use `###`.
- Do not write a `参考文献` section by hand; it is generated from `:cite[...]`.
- The chapter ends with 演習問題 (exercises) and 確認テスト (quiz).

## 3. Language

- である調 (plain form). No casual tone.
- Use standard Japanese mathematical/physical terminology (文部省『学術用語集』準拠). Names of people stay in Latin script: Cauchy, Hilbert, Schrödinger.
- On first definition, introduce terms with `:term[固有値]{en="eigenvalue" yomi="こゆうち"}`. `yomi` must be kana. This feeds the glossary.
- Full-width punctuation 「、」「。」 in prose; ASCII inside math.

## 4. Rigor

- Every theorem/proposition/lemma/corollary has a proof, or an explicit pointer: 「証明は :cite[key] を参照」 with a reason for omission.
- Definitions precede use. Notation is introduced before it is used.
- State hypotheses precisely (field K, finite dimensionality, regularity classes, etc.).
- Physics: state assumptions and approximations explicitly; derivations must be complete enough to reproduce.

## 5. Markup reference

### Environments (numbered `章.通し番号`, shared counter)

```markdown
:::definition[固有値]{#def-linalg-eigenvalue}
...
:::

:::theorem[スペクトル定理]{#thm-linalg-spectral}
...
:::

:::proof
...
:::
```

Kinds and label prefixes:

| directive | 表示 | label prefix |
|---|---|---|
| definition | 定義 | def |
| theorem | 定理 | thm |
| proposition | 命題 | prop |
| lemma | 補題 | lem |
| corollary | 系 | cor |
| example | 例 | exm |
| remark | 注意 (unnumbered) | rem |
| axiom | 公理 | ax |
| postulate | 公準 | post |
| law | 法則 | law |
| principle | 原理 | prin |

Labels are `<prefix>-<subject prefix>-<name>` in lowercase kebab-case, e.g. `thm-linalg-spectral`. Labels are optional; add them when the item is referenced.

`:::proof[定理 5.3 の証明]` customizes the proof heading (use when the proof is separated from its statement).

Nesting: an outer container needs more colons than the inner one (`::::exercise` containing `:::solution`).

### Cross references

- `:ref[thm-linalg-spectral]` → 「定理 5.3」 (prefixed by subject title if in another subject).
- `:ref[thm-linalg-spectral]{full}` → 「定理 5.3（スペクトル定理）」.
- `:ref[thm-linalg-spectral]{text="この定理"}` → custom link text.
- Plain links to subject/chapter pages: `[線形代数](/math/linear-algebra/)` (absolute paths without base).

### Equations

Number an equation by adding `\label{...}` inside display math; reference with `:eqref[...]`.

```markdown
$$
\Phi_A(t) = \det(tI - A) \label{eq-linalg-charpoly}
$$
```

One `\label` per display. Use `aligned` for multi-line equations. Unlabeled equations are unnumbered.

### Exercises

```markdown
::::exercise{level=basic}
問題文
:::hint
ヒント（任意）
:::
:::solution
解答
:::
::::
```

`level` is `basic` (基本) or `advanced` (発展). Solutions must be complete, not just final answers.

### Quiz (確認テスト)

```markdown
::::quiz
問題文

- [ ] 選択肢
- [x] 正解の選択肢
- [ ] 選択肢

:::explanation
解説
:::
::::

::::quiz{answer="3" tol="1e-6"}
数値で答える問題
:::explanation
解説
:::
::::
```

More than one `[x]` makes it a multiple-choice question. Options are labeled (A), (B), … automatically. Every quiz renders a collapsible 模範解答 block (correct option letters or numeric answer, followed by the explanation) that readers can open without answering, so the explanation must be a complete model answer, not just a hint.

### Citations

`:cite[saito-linalg]` or `:cite[saito-linalg]{loc="第4章"}`. Keys are defined in `curriculum/references.yaml`. Add only well-known textbooks; never guess years, editions or page numbers.

### Figures

- Diagrams: fenced code with `tikz` or `tikzcd` (rendered to SVG at build time). Optional caption: ` ```tikzcd caption="..." `.
- Interactive figures: `::interactive[キャプション]{name="pendulum-phase"}` with `src/interactive/<name>.ts` exporting `default (el: HTMLElement, params) => void`.

## 6. Notation and conventions

KaTeX macros (`src/lib/katex-macros.ts`):

| macro | output |
|---|---|
| `\R \C \Z \Q \N \K \F` | blackboard bold sets |
| `\abs{x}` `\norm{x}` `\inner{x}{y}` | \|x\|, ‖x‖, ⟨x, y⟩ |
| `\bra{a}` `\ket{a}` `\braket{a}{b}` `\ketbra{a}{b}` `\mel{a}{A}{b}` `\expval{A}` | Dirac notation |
| `\comm{A}{B}` `\acomm{A}{B}` | [A, B], {A, B} |
| `\dd` `\dv{f}{x}` `\dvn{f}{x}{n}` `\pdv{f}{x}` `\pdvn{f}{x}{n}` `\fdv{F}{\phi}` | derivatives |
| `\vb{r}` | bold vector |
| `\grad \divg \rot` | grad, div, rot |
| `\tr \Tr \rank \diag \sgn \spn \Ker \Img \Re \Im \Hom \End \Aut` | operators |
| `\GL \SL \SO \SU \id` | groups, identity |
| `\eps \defeq \order{x} \conj{z} \transp{A} \adj{A}` | ε, ≔, O(x), z̄, Aᵀ, A* |

Conventions declared per subject (see each subject page):

- 線形代数・関数解析: complex inner product linear in the **first** argument.
- 量子力学: Dirac notation, linear in the **second** argument (ket). Mention the difference when relevant.
- 電磁気学: SI units.
- 特殊相対論・一般相対論: metric signature (−,+,+,+).
- 場の量子論・標準模型: ℏ = c = 1, metric signature (+,−,−,−).

## 7. Self-review checklist

- [ ] All registry sections are covered in order.
- [ ] Definitions precede use; every theorem has a proof or an explicit citation.
- [ ] Non-trivial computations verified (sympy or by hand twice).
- [ ] At least 4 exercises (basic and advanced) with complete solutions, solvable using the chapter and its prerequisites.
- [ ] 3–6 quiz questions; answers and explanations agree.
- [ ] Terms introduced with `:term`; key results labeled; notation listed in frontmatter.
- [ ] `npm run check` and `npm run build` pass.
