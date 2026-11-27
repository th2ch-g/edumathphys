import { readdirSync, readFileSync } from 'node:fs';
import { join, relative, sep } from 'node:path';
import { visit } from 'unist-util-visit';
import { toString } from 'mdast-util-to-string';
import { CONTENT_DIR, loadCurriculum, loadReferences, type Chapter } from './curriculum.ts';
import {
  ENV_KINDS,
  EQUATION_PREFIX,
  EXERCISE_PREFIX,
  numberTree,
  parseMarkdown,
  termAnchor,
  type EnvKind,
} from './markdown.ts';

export interface LabelEntry {
  id: string;
  kind: EnvKind | 'exercise' | 'equation';
  kindLabel: string;
  number?: string;
  title?: string;
  subjectId?: string;
  subjectTitle?: string;
  url: string;
  file: string;
}

export interface TermEntry {
  term: string;
  en: string;
  yomi: string;
  url: string;
  docTitle: string;
  subjectTitle?: string;
}

export interface NotationEntry {
  symbol: string;
  meaning: string;
  url: string;
  docTitle: string;
  subjectId?: string;
}

export interface DocInfo {
  file: string;
  url: string;
  title: string;
  chapter?: Chapter;
  quizIds: string[];
  exerciseIds: string[];
}

export interface DocIndex {
  labels: Map<string, LabelEntry>;
  terms: TermEntry[];
  notation: NotationEntry[];
  docs: Map<string, DocInfo>;
  errors: string[];
}

function listMarkdown(dir: string): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((d) => {
    const p = join(dir, d.name);
    if (d.isDirectory()) return listMarkdown(p);
    return d.name.endsWith('.md') ? [p] : [];
  });
}

export function docUrl(file: string): string {
  const rel = relative(CONTENT_DIR, file).split(sep).join('/').replace(/\.md$/, '');
  if (rel === 'index') return '/';
  return '/' + rel.replace(/\/index$/, '') + '/';
}

/** Progress-tracking key for a document (stable across renumbering). */
export function docKey(file: string): string {
  return docUrl(file).replace(/^\/|\/$/g, '') || 'index';
}

const LABEL_ID_RE = /^[a-z]+-[a-z][a-z0-9]*-[a-z0-9]+(-[a-z0-9]+)*$/;
const KANA_RE = /^[\p{Script=Hiragana}\p{Script=Katakana}ー・\s]+$/u;

let cached: DocIndex | undefined;

export function buildDocIndex(opts: { fresh?: boolean } = {}): DocIndex {
  if (cached && !opts.fresh) return cached;
  const cur = loadCurriculum({ fresh: opts.fresh });
  const { refs } = loadReferences();
  const labels = new Map<string, LabelEntry>();
  const terms: TermEntry[] = [];
  const notation: NotationEntry[] = [];
  const docs = new Map<string, DocInfo>();
  const errors: string[] = [];
  const trees = new Map<string, ReturnType<typeof parseMarkdown>['tree']>();

  for (const file of listMarkdown(CONTENT_DIR)) {
    const rel = relative(process.cwd(), file);
    const chapter = cur.chapterByFile.get(file);
    const subject = chapter?.subject;
    const { tree, frontmatter } = parseMarkdown(readFileSync(file, 'utf8'));
    trees.set(file, tree);
    const url = docUrl(file);
    const title = String(frontmatter.title ?? '');
    const key = docKey(file);
    const info: DocInfo = { file, url, title, chapter, quizIds: [], exerciseIds: [] };
    docs.set(file, info);

    const relParts = relative(CONTENT_DIR, file).split(sep);
    if (!chapter && (relParts[0] === 'math' || relParts[0] === 'physics')) {
      errors.push(`${rel}: not registered in curriculum`);
    }
    if (chapter && title !== chapter.title) {
      errors.push(`${rel}: frontmatter title "${title}" differs from curriculum "${chapter.title}"`);
    }

    const add = (entry: LabelEntry, expectedPrefix: string) => {
      if (!LABEL_ID_RE.test(entry.id) || !entry.id.startsWith(`${expectedPrefix}-`)) {
        errors.push(`${rel}: label "${entry.id}" must look like "${expectedPrefix}-<subject prefix>-<name>" (lowercase, hyphen-separated)`);
      } else if (subject && !entry.id.startsWith(`${expectedPrefix}-${subject.prefix}-`)) {
        errors.push(`${rel}: label "${entry.id}" must be namespaced as "${expectedPrefix}-${subject.prefix}-<name>"`);
      }
      if (labels.has(entry.id)) {
        errors.push(`${rel}: duplicate label "${entry.id}" (also in ${relative(process.cwd(), labels.get(entry.id)!.file)})`);
      }
      labels.set(entry.id, entry);
    };

    const base = { subjectId: subject?.id, subjectTitle: subject?.title, file };
    for (const item of numberTree(tree, chapter?.number)) {
      if (item.type === 'env' && item.label) {
        add(
          { ...base, id: item.label, kind: item.kind, kindLabel: ENV_KINDS[item.kind].label, number: item.number, title: item.title, url: `${url}#${item.label}` },
          ENV_KINDS[item.kind].prefix,
        );
      } else if (item.type === 'exercise') {
        info.exerciseIds.push(`${key}#exr${item.index}`);
        if (item.label) {
          add({ ...base, id: item.label, kind: 'exercise', kindLabel: '演習', number: item.number, url: `${url}#${item.label}` }, EXERCISE_PREFIX);
        }
      } else if (item.type === 'quiz') {
        info.quizIds.push(`${key}#q${item.index}`);
      } else if (item.type === 'equation') {
        add({ ...base, id: item.label, kind: 'equation', kindLabel: '式', number: item.number, url: `${url}#${item.label}` }, EQUATION_PREFIX);
      }
    }

    const seen = new Map<string, number>();
    visit(tree, 'textDirective', (node) => {
      const attrs = node.attributes ?? {};
      if (node.name === 'term') {
        const term = toString(node);
        const en = attrs.en ?? '';
        const yomi = attrs.yomi ?? '';
        if (!term || !en || !yomi) errors.push(`${rel}: :term[${term}] requires en and yomi attributes`);
        else if (!KANA_RE.test(yomi)) errors.push(`${rel}: :term[${term}] yomi "${yomi}" must be kana`);
        const anchor = termAnchor(en || term, seen);
        terms.push({ term, en, yomi, url: `${url}#${anchor}`, docTitle: title, subjectTitle: subject?.title });
      } else if (node.name === 'cite') {
        const k = toString(node).trim();
        if (!refs.has(k)) errors.push(`${rel}: unknown citation key "${k}"`);
      }
    });

    const notes = frontmatter.notation;
    if (notes !== undefined) {
      if (!Array.isArray(notes)) errors.push(`${rel}: frontmatter notation must be a list`);
      else {
        for (const n of notes as { symbol?: unknown; meaning?: unknown }[]) {
          if (typeof n?.symbol !== 'string' || typeof n?.meaning !== 'string') {
            errors.push(`${rel}: notation entries need string symbol and meaning`);
            continue;
          }
          notation.push({ symbol: n.symbol, meaning: n.meaning, url, docTitle: title, subjectId: subject?.id });
        }
      }
    }
  }

  // Validate references after all labels are known.
  for (const [file, tree] of trees) {
    const rel = relative(process.cwd(), file);
    visit(tree, 'textDirective', (node) => {
      if (node.name !== 'ref' && node.name !== 'eqref') return;
      const id = toString(node).trim();
      const entry = labels.get(id);
      if (!entry) errors.push(`${rel}: unresolved :${node.name}[${id}]`);
      else if (node.name === 'eqref' && entry.kind !== 'equation') errors.push(`${rel}: :eqref must point to an equation (${id})`);
    });
  }

  cached = { labels, terms, notation, docs, errors };
  return cached;
}

export function invalidateDocIndex(): void {
  cached = undefined;
}
