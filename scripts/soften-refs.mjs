/**
 * Replace :ref/:eqref directives that point to unpublished chapters with plain links,
 * so that chapters can be published before the chapters they reference.
 *
 * Usage: node scripts/soften-refs.mjs <chapter.md> [...]
 * "Published" = tracked in git HEAD, or one of the files given on the command line.
 */
import { execFileSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import { relative, resolve } from 'node:path';

const root = process.cwd();
const files = process.argv.slice(2).map((f) => resolve(root, f));
const { buildDocIndex } = await import(resolve(root, 'src/lib/doc-index.ts'));

const tracked = new Set(
  execFileSync('git', ['ls-tree', '-r', '--name-only', 'HEAD', 'src/content/docs'], { cwd: root, encoding: 'utf8' })
    .split('\n')
    .filter(Boolean)
    .map((f) => resolve(root, f)),
);
const published = new Set([...tracked, ...files]);
const { labels } = buildDocIndex();

for (const file of files) {
  let text = readFileSync(file, 'utf8');
  let count = 0;
  text = text.replace(/:(ref|eqref)\[([a-z0-9-]+)\](\{[^}]*\})?/g, (m, kind, id) => {
    const entry = labels.get(id);
    if (!entry || published.has(entry.file)) return m;
    count++;
    const what = entry.title ?? (kind === 'eqref' ? '式' : entry.kindLabel);
    const where = entry.subjectTitle ? `${entry.subjectTitle}の` : '';
    return `[${where}${what}](${entry.url})`;
  });
  if (count > 0) {
    writeFileSync(file, text);
    console.log(`${relative(root, file)}: softened ${count} reference(s)`);
  }
}
