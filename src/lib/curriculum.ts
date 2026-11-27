import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { parse as parseYaml } from 'yaml';
import { z } from 'zod';

export const ROOT = process.cwd();
export const CURRICULUM_DIR = join(ROOT, 'curriculum');
export const CONTENT_DIR = join(ROOT, 'src', 'content', 'docs');

export const AREAS = { math: '数学', physics: '物理' } as const;
export type Area = keyof typeof AREAS;

export const LEVELS = {
  'undergrad-basic': '学部基礎',
  'undergrad-advanced': '学部専門',
  graduate: '大学院',
} as const;
export type Level = keyof typeof LEVELS;

export const STATUSES = { stub: '執筆中', draft: '草稿', reviewed: '公開' } as const;
export type Status = keyof typeof STATUSES;

const slug = z.string().regex(/^[a-z0-9]+(-[a-z0-9]+)*$/, 'must be kebab-case');

const chapterSchema = z
  .object({
    id: slug,
    title: z.string().min(1),
    status: z.enum(['stub', 'draft', 'reviewed']).default('stub'),
    prerequisites: z.array(z.string()).default([]),
    sections: z.array(z.string().min(1)).min(1),
  })
  .strict();

const subjectSchema = z
  .object({
    id: slug,
    prefix: z.string().regex(/^[a-z][a-z0-9]*$/),
    title: z.string().min(1),
    title_en: z.string().min(1),
    area: z.enum(['math', 'physics']),
    level: z.enum(['undergrad-basic', 'undergrad-advanced', 'graduate']),
    order: z.number().int(),
    description: z.string().min(1),
    prerequisites: z.array(slug).default([]),
    conventions: z.array(z.string()).default([]),
    references: z.array(z.string()).default([]),
    chapters: z.array(chapterSchema).min(1),
  })
  .strict();

export type ChapterDef = z.infer<typeof chapterSchema>;
export type SubjectDef = z.infer<typeof subjectSchema>;

export interface Chapter extends ChapterDef {
  key: string;
  number: number;
  subject: Subject;
  url: string;
  file: string;
  exists: boolean;
}

export interface Subject extends Omit<SubjectDef, 'chapters'> {
  chapters: Chapter[];
  url: string;
  sourceFile: string;
}

export interface Curriculum {
  subjects: Subject[];
  subjectById: Map<string, Subject>;
  subjectByPrefix: Map<string, Subject>;
  chapterByKey: Map<string, Chapter>;
  chapterByFile: Map<string, Chapter>;
  errors: string[];
}

export function withBase(path: string): string {
  const base = (process.env.SITE_BASE ?? '/').replace(/\/?$/, '/');
  return base === '/' ? path : base + path.replace(/^\//, '');
}

function listYaml(dir: string): string[] {
  if (!existsSync(dir)) return [];
  return readdirSync(dir, { withFileTypes: true }).flatMap((d) => {
    const p = join(dir, d.name);
    if (d.isDirectory()) return listYaml(p);
    return d.name.endsWith('.yaml') ? [p] : [];
  });
}

let cached: Curriculum | undefined;

export function loadCurriculum(opts: { fresh?: boolean } = {}): Curriculum {
  if (cached && !opts.fresh) return cached;
  const errors: string[] = [];
  const subjects: Subject[] = [];
  for (const file of listYaml(join(CURRICULUM_DIR, 'subjects'))) {
    const rel = file.slice(ROOT.length + 1);
    const parsed = subjectSchema.safeParse(parseYaml(readFileSync(file, 'utf8')));
    if (!parsed.success) {
      errors.push(`${rel}: ${parsed.error.issues.map((i) => `${i.path.join('.')}: ${i.message}`).join('; ')}`);
      continue;
    }
    const def = parsed.data;
    if (!rel.endsWith(`${def.area}/${def.id}.yaml`)) {
      errors.push(`${rel}: file must be located at curriculum/subjects/${def.area}/${def.id}.yaml`);
    }
    const subject: Subject = {
      ...def,
      chapters: [],
      url: `/${def.area}/${def.id}/`,
      sourceFile: rel,
    };
    subject.chapters = def.chapters.map((c, i) => {
      const filePath = join(CONTENT_DIR, def.area, def.id, `${c.id}.md`);
      return {
        ...c,
        key: `${def.id}/${c.id}`,
        number: i + 1,
        subject,
        url: `/${def.area}/${def.id}/${c.id}/`,
        file: filePath,
        exists: existsSync(filePath),
      };
    });
    subjects.push(subject);
  }
  subjects.sort((a, b) => (a.area === b.area ? a.order - b.order : a.area === 'math' ? -1 : 1));

  const subjectById = new Map<string, Subject>();
  const subjectByPrefix = new Map<string, Subject>();
  const chapterByKey = new Map<string, Chapter>();
  const chapterByFile = new Map<string, Chapter>();
  for (const s of subjects) {
    if (subjectById.has(s.id)) errors.push(`duplicate subject id: ${s.id}`);
    if (subjectByPrefix.has(s.prefix)) errors.push(`duplicate subject prefix: ${s.prefix} (${s.id})`);
    subjectById.set(s.id, s);
    subjectByPrefix.set(s.prefix, s);
    for (const c of s.chapters) {
      if (chapterByKey.has(c.key)) errors.push(`duplicate chapter id: ${c.key}`);
      chapterByKey.set(c.key, c);
      chapterByFile.set(c.file, c);
      if (c.exists && c.status === 'stub') errors.push(`${c.key}: file exists but status is stub`);
      if (!c.exists && c.status !== 'stub') errors.push(`${c.key}: status is ${c.status} but ${c.file.slice(ROOT.length + 1)} is missing`);
    }
  }

  for (const s of subjects) {
    for (const p of s.prerequisites) {
      if (!subjectById.has(p)) errors.push(`${s.id}: unknown prerequisite subject "${p}"`);
    }
    for (const c of s.chapters) {
      for (const p of c.prerequisites) {
        if (!chapterByKey.has(p) && !subjectById.has(p)) errors.push(`${c.key}: unknown prerequisite "${p}"`);
      }
    }
  }
  const cycle = findCycle(subjects, subjectById);
  if (cycle) errors.push(`prerequisite cycle: ${cycle.join(' -> ')}`);

  cached = { subjects, subjectById, subjectByPrefix, chapterByKey, chapterByFile, errors };
  return cached;
}

/** Subject-level dependency edges, including those implied by chapter-level prerequisites. */
export function subjectDependencies(s: Subject, byId: Map<string, Subject>): Set<string> {
  const deps = new Set(s.prerequisites);
  for (const c of s.chapters) {
    for (const p of c.prerequisites) {
      const dep = p.split('/')[0];
      if (dep !== s.id && byId.has(dep)) deps.add(dep);
    }
  }
  return deps;
}

function findCycle(subjects: Subject[], byId: Map<string, Subject>): string[] | undefined {
  const state = new Map<string, 'visiting' | 'done'>();
  const stack: string[] = [];
  const visit = (id: string): string[] | undefined => {
    const st = state.get(id);
    if (st === 'done') return undefined;
    if (st === 'visiting') return [...stack.slice(stack.indexOf(id)), id];
    state.set(id, 'visiting');
    stack.push(id);
    const s = byId.get(id);
    if (s) {
      for (const d of subjectDependencies(s, byId)) {
        const c = visit(d);
        if (c) return c;
      }
    }
    stack.pop();
    state.set(id, 'done');
    return undefined;
  };
  for (const s of subjects) {
    const c = visit(s.id);
    if (c) return c;
  }
  return undefined;
}

export function prerequisiteLinks(items: string[], cur: Curriculum): { label: string; url: string }[] {
  return items.flatMap((p) => {
    const ch = cur.chapterByKey.get(p);
    if (ch) return [{ label: `${ch.subject.title} 第${ch.number}章 ${ch.title}`, url: ch.url }];
    const s = cur.subjectById.get(p);
    return s ? [{ label: s.title, url: s.url }] : [];
  });
}

const referenceSchema = z
  .object({
    author: z.string().min(1),
    title: z.string().min(1),
    publisher: z.string().min(1),
    year: z.number().int().optional(),
    edition: z.string().optional(),
    lang: z.enum(['ja', 'en']),
  })
  .strict();
export type Reference = z.infer<typeof referenceSchema> & { key: string };

let cachedRefs: { refs: Map<string, Reference>; errors: string[] } | undefined;

export function loadReferences(opts: { fresh?: boolean } = {}): { refs: Map<string, Reference>; errors: string[] } {
  if (cachedRefs && !opts.fresh) return cachedRefs;
  const errors: string[] = [];
  const refs = new Map<string, Reference>();
  const file = join(CURRICULUM_DIR, 'references.yaml');
  const raw = existsSync(file) ? (parseYaml(readFileSync(file, 'utf8')) ?? {}) : {};
  for (const [key, value] of Object.entries(raw as Record<string, unknown>)) {
    const parsed = referenceSchema.safeParse(value);
    if (parsed.success) refs.set(key, { ...parsed.data, key });
    else errors.push(`references.yaml ${key}: ${parsed.error.issues.map((i) => `${i.path.join('.')}: ${i.message}`).join('; ')}`);
  }
  cachedRefs = { refs, errors };
  return cachedRefs;
}

export function formatReference(r: Reference): string {
  const ed = r.edition ? `${r.edition}, ` : '';
  const year = r.year ? `, ${r.year}` : '';
  return r.lang === 'ja'
    ? `${r.author}『${r.title}』${ed}${r.publisher}${year}.`
    : `${r.author}, *${r.title}*, ${ed}${r.publisher}${year}.`;
}
