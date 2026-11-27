import type { Root, Nodes, Parent, Paragraph } from 'mdast';
import type { ContainerDirective } from 'mdast-util-directive';
import { toString } from 'mdast-util-to-string';
import remarkDirective from 'remark-directive';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import remarkParse from 'remark-parse';
import { parse as parseYaml } from 'yaml';
import { unified } from 'unified';
import { visit } from 'unist-util-visit';

export const ENV_KINDS = {
  definition: { label: '定義', prefix: 'def' },
  theorem: { label: '定理', prefix: 'thm' },
  proposition: { label: '命題', prefix: 'prop' },
  lemma: { label: '補題', prefix: 'lem' },
  corollary: { label: '系', prefix: 'cor' },
  example: { label: '例', prefix: 'exm' },
  remark: { label: '注意', prefix: 'rem' },
  axiom: { label: '公理', prefix: 'ax' },
  postulate: { label: '公準', prefix: 'post' },
  law: { label: '法則', prefix: 'law' },
  principle: { label: '原理', prefix: 'prin' },
} as const;
export type EnvKind = keyof typeof ENV_KINDS;
export const UNNUMBERED_ENVS = new Set<EnvKind>(['remark']);

export const EXERCISE_PREFIX = 'exr';
export const EQUATION_PREFIX = 'eq';

export const BLOCK_DIRECTIVES = new Set([
  ...Object.keys(ENV_KINDS),
  'proof',
  'exercise',
  'solution',
  'hint',
  'quiz',
  'explanation',
  'interactive',
  // Starlight asides
  'note',
  'tip',
  'caution',
  'danger',
]);
export const TEXT_DIRECTIVES = new Set(['ref', 'eqref', 'cite', 'term']);

export const EXERCISE_LEVELS = { basic: '基本', advanced: '発展' } as const;

const LABEL_RE = /\\label\{([^}]+)\}/g;

export function isEnvKind(name: string): name is EnvKind {
  return Object.hasOwn(ENV_KINDS, name);
}

export function directiveLabel(node: ContainerDirective): Paragraph | undefined {
  const first = node.children[0];
  if (first?.type === 'paragraph' && (first.data as { directiveLabel?: boolean } | undefined)?.directiveLabel) return first;
  return undefined;
}

export function extractEquationLabel(value: string): { label?: string; value: string; error?: string } {
  const labels = [...value.matchAll(LABEL_RE)].map((m) => m[1].trim());
  if (labels.length === 0) return { value };
  const stripped = value.replace(LABEL_RE, '');
  if (labels.length > 1) return { value: stripped, label: labels[0], error: `multiple \\label in one display math: ${labels.join(', ')}` };
  return { label: labels[0], value: stripped };
}

export type NumberedItem =
  | { type: 'env'; node: ContainerDirective; kind: EnvKind; number?: string; label?: string; title?: string }
  | { type: 'exercise'; node: ContainerDirective; number: string; index: number; label?: string }
  | { type: 'quiz'; node: ContainerDirective; number: string; index: number }
  | { type: 'equation'; node: Nodes & { type: 'math' }; number: string; label: string };

/**
 * Assigns numbers in document order. Both the index scanner and the renderer call this,
 * so numbering is guaranteed to be identical.
 */
export function numberTree(tree: Root, chapterNumber?: number): NumberedItem[] {
  const pre = chapterNumber ? `${chapterNumber}.` : '';
  let env = 0;
  let exercise = 0;
  let quiz = 0;
  let equation = 0;
  const items: NumberedItem[] = [];
  visit(tree, (node) => {
    if (node.type === 'containerDirective') {
      const name = node.name;
      const id = typeof node.attributes?.id === 'string' ? node.attributes.id : undefined;
      if (isEnvKind(name)) {
        const labelPara = directiveLabel(node);
        const title = labelPara ? toString(labelPara) : undefined;
        if (UNNUMBERED_ENVS.has(name)) items.push({ type: 'env', node, kind: name, label: id, title });
        else items.push({ type: 'env', node, kind: name, number: `${pre}${++env}`, label: id, title });
      } else if (name === 'exercise') {
        exercise++;
        items.push({ type: 'exercise', node, number: `${pre}${exercise}`, index: exercise, label: id });
      } else if (name === 'quiz') {
        quiz++;
        items.push({ type: 'quiz', node, number: `${quiz}`, index: quiz });
      }
    } else if (node.type === 'math') {
      const { label } = extractEquationLabel(node.value);
      if (label) items.push({ type: 'equation', node, number: `${pre}${++equation}`, label });
    }
  });
  return items;
}

export function parseMarkdown(source: string): { tree: Root; frontmatter: Record<string, unknown> } {
  let body = source;
  let frontmatter: Record<string, unknown> = {};
  const m = /^---\r?\n([\s\S]*?)\r?\n---\r?\n?/.exec(source);
  if (m) {
    frontmatter = (parseYaml(m[1]) as Record<string, unknown>) ?? {};
    body = source.slice(m[0].length);
  }
  const processor = unified().use(remarkParse).use(remarkGfm).use(remarkMath).use(remarkDirective);
  const tree = processor.runSync(processor.parse(body)) as Root;
  return { tree, frontmatter };
}

export function termAnchor(en: string, seen: Map<string, number>): string {
  const base = 'term-' + (en.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '') || 'x');
  const n = (seen.get(base) ?? 0) + 1;
  seen.set(base, n);
  return n === 1 ? base : `${base}-${n}`;
}

export function isParent(node: Nodes): node is Nodes & Parent {
  return 'children' in node;
}
