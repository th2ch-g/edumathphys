import { existsSync } from 'node:fs';
import type {} from 'mdast-util-to-hast';
import { join } from 'node:path';
import type { Root, Nodes, Parent, PhrasingContent, RootContent, List, ListItem, Paragraph, Code } from 'mdast';
import type { ContainerDirective, LeafDirective, TextDirective } from 'mdast-util-directive';
import { toString } from 'mdast-util-to-string';
import { visit, SKIP } from 'unist-util-visit';
import type { VFile } from 'vfile';
import { loadCurriculum, loadReferences, withBase, type Reference } from '../lib/curriculum.ts';
import { buildDocIndex, docKey, type LabelEntry } from '../lib/doc-index.ts';
import {
  BLOCK_DIRECTIVES,
  ENV_KINDS,
  EXERCISE_LEVELS,
  TEXT_DIRECTIVES,
  directiveLabel,
  extractEquationLabel,
  numberTree,
  termAnchor,
} from '../lib/markdown.ts';
import { renderTikz } from './tikz.ts';

const INTERACTIVE_DIR = join(process.cwd(), 'src', 'interactive');

type HData = { hName?: string; hProperties?: Record<string, unknown>; hChildren?: unknown[] };

function setH(node: Nodes, hName: string, hProperties: Record<string, unknown> = {}): void {
  const data = (node.data ??= {}) as HData;
  data.hName = hName;
  data.hProperties = { ...(data.hProperties ?? {}), ...hProperties };
}

function span(className: string, children: PhrasingContent[]): PhrasingContent {
  return { type: 'emphasis', data: { hName: 'span', hProperties: { className: [className] } }, children } as PhrasingContent;
}

/** remark-math stores the TeX source both in `value` and in prebuilt `hChildren`; keep them in sync. */
function setMathValue(node: { value: string; data?: unknown }, value: string): void {
  node.value = value;
  const code = (node.data as HData | undefined)?.hChildren?.[0] as { children?: { type: string; value: string }[] } | undefined;
  const leaf = code?.children?.[0];
  if (leaf?.type === 'text') leaf.value = value;
}

function text(value: string): PhrasingContent {
  return { type: 'text', value };
}

function headParagraph(className: string, children: PhrasingContent[]): Paragraph {
  return { type: 'paragraph', data: { hName: 'p', hProperties: { className: [className] } }, children };
}

function refText(entry: LabelEntry, currentSubject: string | undefined): string {
  const prefix = entry.subjectId && entry.subjectId !== currentSubject && entry.subjectTitle ? `${entry.subjectTitle} ` : '';
  if (entry.kind === 'equation') return `${prefix}(${entry.number})`;
  return `${prefix}${entry.kindLabel}${entry.number ? ` ${entry.number}` : ''}`;
}

function referenceHtmlChildren(r: Reference): PhrasingContent[] {
  const ed = r.edition ? `${r.edition}, ` : '';
  const year = r.year ? `, ${r.year}` : '';
  if (r.lang === 'ja') return [text(`${r.author}『${r.title}』${ed}${r.publisher}${year}.`)];
  return [text(`${r.author}, `), { type: 'emphasis', children: [text(r.title)] }, text(`, ${ed}${r.publisher}${year}.`)];
}

export interface RemarkEduMathOptions {
  /** Only files under these absolute directory prefixes are transformed. */
  include: string[];
}

export default function remarkEduMath(options: RemarkEduMathOptions) {
  return async function transformer(tree: Root, file: VFile): Promise<void> {
    const path = file.path;
    if (!path || !options.include.some((p) => path.startsWith(p))) return;
    const cur = loadCurriculum();
    const index = buildDocIndex();
    const { refs } = loadReferences();
    const chapter = cur.chapterByFile.get(path);
    const subjectId = chapter?.subject.id;
    const key = docKey(path);
    const fail = (msg: string, node?: Nodes): never => file.fail(msg, node ? { place: node.position } : undefined);
    const labeledEquations = new Map<Nodes, string>();

    // 1. Unknown directives are almost always typos; fail loudly.
    visit(tree, (node) => {
      if ((node.type === 'containerDirective' || node.type === 'leafDirective') && !BLOCK_DIRECTIVES.has(node.name)) {
        fail(`unknown directive "${node.name}"`, node);
      }
      if (node.type === 'textDirective' && !TEXT_DIRECTIVES.has(node.name) && (node.children.length > 0 || Object.keys(node.attributes ?? {}).length > 0)) {
        fail(`unknown inline directive ":${node.name}"`, node);
      }
    });

    // 2. Numbered blocks (numbering shared with the index scanner).
    for (const item of numberTree(tree, chapter?.number)) {
      if (item.type === 'env') {
        const node = item.node;
        const labelPara = directiveLabel(node);
        const headChildren: PhrasingContent[] = [span('env-label', [text(`${ENV_KINDS[item.kind].label}${item.number ? ` ${item.number}` : ''}`)])];
        if (labelPara) headChildren.push(span('env-title', [text('（'), ...labelPara.children, text('）')]));
        const head = headParagraph('env-head', headChildren);
        if (labelPara) node.children.splice(0, 1, head);
        else node.children.unshift(head);
        setH(node, 'div', { className: ['env', `env-${item.kind}`], ...(item.label ? { id: item.label } : {}) });
      } else if (item.type === 'exercise') {
        const node = item.node;
        const level = (node.attributes?.level ?? 'basic') as keyof typeof EXERCISE_LEVELS;
        if (!Object.hasOwn(EXERCISE_LEVELS, level)) fail(`exercise level must be basic or advanced (got "${level}")`, node);
        const labelPara = directiveLabel(node);
        const headChildren: PhrasingContent[] = [span('env-label', [text(`演習 ${item.number}`)]), span(`badge badge-${level}`, [text(EXERCISE_LEVELS[level])])];
        if (labelPara) headChildren.push(span('env-title', [text('（'), ...labelPara.children, text('）')]));
        const head = headParagraph('env-head', headChildren);
        if (labelPara) node.children.splice(0, 1, head);
        else node.children.unshift(head);
        setH(node, 'div', { className: ['exercise'], id: item.label ?? `exr-${item.index}`, 'data-exid': `${key}#exr${item.index}` });
      } else if (item.type === 'quiz') {
        transformQuiz(item.node, item.number, `${key}#q${item.index}`, fail);
      } else if (item.type === 'equation') {
        const node = item.node as Nodes & { type: 'math'; value: string };
        const { value, error } = extractEquationLabel(node.value);
        if (error) fail(error, node);
        setMathValue(node, `${value.trim()}\n\\tag{${item.number}}`);
        labeledEquations.set(node, item.label);
      }
    }

    // 3. Remaining block directives and equation wrappers.
    visit(tree, (node, _idx, parent) => {
      if (node.type !== 'containerDirective' && node.type !== 'leafDirective') return;
      if (node.type === 'containerDirective' && (node.name === 'proof')) {
        transformProof(node);
      } else if (node.type === 'containerDirective' && (node.name === 'solution' || node.name === 'hint')) {
        if (parent?.type !== 'containerDirective' || (parent as ContainerDirective).name !== 'exercise') {
          if (node.name === 'solution') fail(':::solution must be nested inside ::::exercise', node);
        }
        const labelPara = directiveLabel(node);
        const summary: Paragraph = { type: 'paragraph', data: { hName: 'summary' }, children: labelPara ? labelPara.children : [text(node.name === 'solution' ? '解答' : 'ヒント')] };
        if (labelPara) node.children.splice(0, 1, summary);
        else node.children.unshift(summary);
        setH(node, 'details', { className: [node.name] });
      } else if (node.type === 'containerDirective' && node.name === 'explanation') {
        if (parent?.type !== 'containerDirective' || (parent as ContainerDirective).name !== 'quiz') fail(':::explanation must be nested inside ::::quiz', node);
        setH(node, 'div', { className: ['quiz-explanation'], hidden: true });
        node.children.unshift(headParagraph('quiz-explanation-head', [text('解説')]));
      } else if (node.type === 'leafDirective' && node.name === 'interactive') {
        transformInteractive(node, fail);
      }
    });

    // Wrap labeled equations so they carry an anchor id.
    visit(tree, 'math', (node, idx, parent) => {
      if (!parent || idx === undefined) return;
      const item = labeledEquations.get(node);
      if (!item) return;
      const wrapper = { type: 'equationWrapper', data: { hName: 'div', hProperties: { id: item, className: ['equation'] } }, children: [node] };
      (parent.children as unknown[]).splice(idx, 1, wrapper);
      return SKIP;
    });

    // 4. Inline directives.
    const cited: string[] = [];
    const seenTerms = new Map<string, number>();
    visit(tree, 'textDirective', (node: TextDirective, idx, parent) => {
      if (!parent || idx === undefined) return;
      const attrs = node.attributes ?? {};
      let replacement: PhrasingContent[] | undefined;
      if (node.name === 'ref' || node.name === 'eqref') {
        const id = toString(node).trim();
        const entry = index.labels.get(id);
        if (!entry) fail(`unresolved :${node.name}[${id}]`, node);
        let label = attrs.text ?? refText(entry!, subjectId);
        if (attrs.full !== undefined && entry!.title) label += `（${entry!.title}）`;
        replacement = [{ type: 'link', url: entry!.url, data: { hProperties: { className: ['xref'], 'data-xref': id } }, children: [text(label)] }];
      } else if (node.name === 'cite') {
        const k = toString(node).trim();
        if (!refs.has(k)) fail(`unknown citation key "${k}"`, node);
        if (!cited.includes(k)) cited.push(k);
        const n = cited.indexOf(k) + 1;
        const loc = attrs.loc ? `, ${attrs.loc}` : '';
        replacement = [{ type: 'link', url: `#bib-${k}`, data: { hProperties: { className: ['cite'] } }, children: [text(`[${n}${loc}]`)] }];
      } else if (node.name === 'term') {
        const term = toString(node);
        const en = attrs.en ?? '';
        if (!en || !attrs.yomi) fail(`:term[${term}] requires en and yomi attributes`, node);
        const anchor = termAnchor(en, seenTerms);
        replacement = [
          { type: 'emphasis', data: { hName: 'dfn', hProperties: { id: anchor, className: ['term'] } }, children: node.children } as PhrasingContent,
          span('term-en', [text(`（${en}）`)]),
        ];
      }
      if (replacement) {
        (parent.children as PhrasingContent[]).splice(idx, 1, ...replacement);
        return [SKIP, idx + replacement.length];
      }
      return undefined;
    });

    if (cited.length > 0) {
      const list: List = {
        type: 'list',
        ordered: true,
        spread: false,
        data: { hProperties: { className: ['bibliography'] } },
        children: cited.map((k): ListItem => ({
          type: 'listItem',
          spread: false,
          data: { hProperties: { id: `bib-${k}` } },
          children: [{ type: 'paragraph', children: referenceHtmlChildren(refs.get(k)!) }],
        })),
      };
      tree.children.push({ type: 'heading', depth: 2, children: [text('参考文献')] }, list);
    }

    // Absolute internal links in content are written without the deploy base path.
    visit(tree, 'link', (node) => {
      if (node.url.startsWith('/') && !node.url.startsWith('//')) node.url = withBase(node.url);
    });

    // 5. TikZ code blocks.
    const tikzNodes: { node: Code; idx: number; parent: Parent }[] = [];
    visit(tree, 'code', (node, idx, parent) => {
      if ((node.lang === 'tikz' || node.lang === 'tikzcd') && parent && idx !== undefined) tikzNodes.push({ node, idx, parent });
    });
    for (const { node, idx, parent } of tikzNodes.reverse()) {
      let svg: string;
      try {
        svg = await renderTikz(node.lang as 'tikz' | 'tikzcd', node.value);
      } catch (e) {
        fail(`TikZ rendering failed: ${(e as Error).message}`, node);
      }
      const caption = /caption="([^"]*)"/.exec(node.meta ?? '')?.[1];
      const html = `<figure class="tikz">${svg!}${caption ? `<figcaption>${escapeHtml(caption)}</figcaption>` : ''}</figure>`;
      parent.children.splice(idx, 1, { type: 'html', value: html } as RootContent);
    }
  };

  function transformProof(node: ContainerDirective): void {
    const labelPara = directiveLabel(node);
    const headText: PhrasingContent[] = labelPara ? labelPara.children : [text('証明')];
    if (labelPara) node.children.shift();
    const head = span('proof-head', [...headText, text('.')]);
    const first = node.children[0];
    if (first?.type === 'paragraph') first.children.unshift(head, text(' '));
    else node.children.unshift(headParagraph('proof-head-block', [head]));
    const last = node.children[node.children.length - 1];
    const qed = span('qed', [text('∎')]);
    if (last?.type === 'paragraph') last.children.push(qed);
    else node.children.push(headParagraph('qed-block', [qed]));
    setH(node, 'div', { className: ['proof'] });
  }

  function transformQuiz(node: ContainerDirective, number: string, qid: string, fail: (m: string, n?: Nodes) => never): void {
    const attrs = node.attributes ?? {};
    const list = node.children.find((c): c is List => c.type === 'list' && c.children.some((li) => typeof li.checked === 'boolean'));
    const props: Record<string, unknown> = { className: ['quiz'], 'data-qid': qid, id: `q${qid.split('#q')[1]}` };
    if (list) {
      const correct = list.children.flatMap((li, i) => (li.checked ? [i] : []));
      if (correct.length === 0) fail('quiz needs at least one correct option ([x])', node);
      if (list.children.some((li) => typeof li.checked !== 'boolean')) fail('every quiz option must be a task item ([ ] or [x])', node);
      const multiple = attrs.multiple !== undefined || correct.length > 1;
      list.children.forEach((li, i) => {
        li.checked = null;
        li.data = { ...(li.data ?? {}), hProperties: { className: ['quiz-option'], 'data-index': String(i) } } as ListItem['data'];
      });
      list.data = { ...(list.data ?? {}), hProperties: { className: ['quiz-options'], role: multiple ? 'group' : 'radiogroup' } } as List['data'];
      Object.assign(props, { 'data-type': multiple ? 'multiple' : 'single', 'data-correct': correct.join(',') });
    } else if (attrs.answer !== undefined) {
      if (Number.isNaN(Number(attrs.answer))) fail(`quiz answer must be numeric (got "${attrs.answer}")`, node);
      Object.assign(props, { 'data-type': 'numeric', 'data-answer': attrs.answer, 'data-tol': attrs.tol ?? '0' });
    } else {
      fail('quiz needs a task list of options or an answer attribute', node);
    }
    node.children.unshift(headParagraph('quiz-head', [span('env-label', [text(`問 ${number}`)])]));
    setH(node, 'div', props);
  }

  function transformInteractive(node: LeafDirective, fail: (m: string, n?: Nodes) => never): void {
    const attrs = node.attributes ?? {};
    const name = attrs.name ?? '';
    if (!/^[a-z0-9-]+$/.test(name) || !existsSync(join(INTERACTIVE_DIR, `${name}.ts`))) {
      fail(`unknown interactive figure "${name}" (expected src/interactive/${name}.ts)`, node);
    }
    const { name: _n, ...params } = attrs;
    const caption = toString(node);
    setH(node, 'figure', { className: ['interactive'], 'data-name': name, 'data-params': JSON.stringify(params) });
    (node.data as HData).hChildren = [
      { type: 'element', tagName: 'div', properties: { className: ['interactive-canvas'] }, children: [{ type: 'element', tagName: 'p', properties: { className: ['interactive-fallback'] }, children: [{ type: 'text', value: 'インタラクティブ図（JavaScript が必要）' }] }] },
      ...(caption ? [{ type: 'element', tagName: 'figcaption', properties: {}, children: [{ type: 'text', value: caption }] }] : []),
    ];
  }
}

function escapeHtml(s: string): string {
  return s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]!);
}
