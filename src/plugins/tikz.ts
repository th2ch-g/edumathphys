import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { join } from 'node:path';

const CACHE_DIR = join(process.cwd(), '.cache', 'tikz');
const RENDER_VERSION = '2';
// TeX points are small next to body text; enlarge so diagram labels match the surrounding font size.
const SCALE = 1.6;

const PACKAGES = { amsmath: '', amssymb: '', 'tikz-cd': '' };
const LIBRARIES = 'arrows.meta,calc,positioning,decorations.pathmorphing,decorations.markings,shapes.geometric';

type Tex2Svg = (source: string, opts: Record<string, unknown>) => Promise<string>;
let tex2svg: Tex2Svg | undefined;
let queue: Promise<unknown> = Promise.resolve();

function load(): Tex2Svg {
  if (!tex2svg) {
    const mod = createRequire(import.meta.url)('node-tikzjax');
    tex2svg = (mod.default ?? mod) as Tex2Svg;
  }
  return tex2svg;
}

export function wrapTikz(lang: 'tikz' | 'tikzcd', code: string): string {
  if (code.includes('\\begin{document}')) return code;
  const env = lang === 'tikzcd' ? 'tikzcd' : 'tikzpicture';
  const body = code.includes(`\\begin{${env}}`) ? code : `\\begin{${env}}\n${code}\n\\end{${env}}`;
  return `\\begin{document}\n${body}\n\\end{document}`;
}

function themeable(svg: string): string {
  return svg
    .replace(/(stroke|fill)="#000(?:000)?"/g, '$1="currentColor"')
    .replace(/<svg ([^>]*)>/, (_m, attrs: string) => {
      const scaled = attrs.replace(/\b(width|height)="([\d.]+)"/g, (_a, k: string, v: string) => `${k}="${(Number(v) * SCALE).toFixed(2)}"`);
      return `<svg fill="currentColor" role="img" ${scaled}>`;
    });
}

/** Renders TikZ to SVG with a content-addressed cache. Calls are serialized (node-tikzjax is not reentrant). */
export function renderTikz(lang: 'tikz' | 'tikzcd', code: string): Promise<string> {
  const source = wrapTikz(lang, code);
  const hash = createHash('sha256').update(RENDER_VERSION).update(LIBRARIES).update(source).digest('hex').slice(0, 24);
  const cacheFile = join(CACHE_DIR, `${hash}.svg`);
  if (existsSync(cacheFile)) return Promise.resolve(readFileSync(cacheFile, 'utf8'));
  const job = queue.then(async () => {
    if (existsSync(cacheFile)) return readFileSync(cacheFile, 'utf8');
    const svg = await load()(source, { texPackages: PACKAGES, tikzLibraries: LIBRARIES });
    if (!svg || !svg.includes('<svg')) throw new Error('TikZ rendering produced no SVG');
    const out = themeable(svg);
    mkdirSync(CACHE_DIR, { recursive: true });
    writeFileSync(cacheFile, out);
    return out;
  });
  queue = job.catch(() => undefined);
  return job;
}
