/**
 * Stage reviewed chapters for a commit while other drafts stay uncommitted.
 *
 * Usage: node scripts/stage-chapters.mjs <subject>/<chapter> [...]
 *
 * 1. Builds the registry as committed (HEAD) with `status: reviewed` set for the given chapters.
 * 2. Verifies `npm run check` and `npm run build` in a temporary worktree containing only
 *    HEAD + the given chapters (so unpublished drafts cannot mask missing dependencies).
 * 3. On success, stages the chapter files and the registry (HEAD + reviewed lines) in the index,
 *    and marks the chapters as reviewed in the working-tree registry too.
 */
import { execFileSync } from 'node:child_process';
import { cpSync, existsSync, mkdtempSync, readFileSync, rmSync, symlinkSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const root = process.cwd();
const keys = process.argv.slice(2);
if (keys.length === 0) {
  console.error('usage: node scripts/stage-chapters.mjs <subject>/<chapter> [...]');
  process.exit(2);
}

const git = (...args) => execFileSync('git', args, { cwd: root, encoding: 'utf8' });

function subjectFile(subject) {
  for (const area of ['math', 'physics']) {
    const p = `curriculum/subjects/${area}/${subject}.yaml`;
    if (existsSync(join(root, p))) return { area, path: p };
  }
  throw new Error(`unknown subject ${subject}`);
}

/** Sets `status: <status>` for chapter `id` in YAML text (keeps the rest byte-identical). */
function setStatus(text, id, status) {
  const lines = text.split('\n');
  const start = lines.findIndex((l) => l.trim() === `- id: ${id}`);
  if (start < 0) throw new Error(`chapter ${id} not found`);
  let end = lines.length;
  for (let i = start + 1; i < lines.length; i++) {
    if (/^\s*- id: /.test(lines[i]) || /^\S/.test(lines[i])) {
      end = i;
      break;
    }
  }
  const statusIdx = lines.slice(start, end).findIndex((l) => /^\s+status:/.test(l));
  if (statusIdx >= 0) lines[start + statusIdx] = `    status: ${status}`;
  else {
    const titleIdx = lines.slice(start, end).findIndex((l) => /^\s+title:/.test(l));
    lines.splice(start + titleIdx + 1, 0, `    status: ${status}`);
  }
  return lines.join('\n');
}

const bySubject = new Map();
for (const key of keys) {
  const [subject, chapter] = key.split('/');
  if (!subject || !chapter) throw new Error(`bad key ${key}`);
  if (!bySubject.has(subject)) bySubject.set(subject, []);
  bySubject.get(subject).push(chapter);
}

const yamlPub = new Map();
const chapterFiles = [];
for (const [subject, chapters] of bySubject) {
  const { area, path } = subjectFile(subject);
  let text;
  try {
    text = git('show', `HEAD:${path}`);
  } catch {
    text = readFileSync(join(root, path), 'utf8');
  }
  for (const c of chapters) {
    text = setStatus(text, c, 'reviewed');
    const f = `src/content/docs/${area}/${subject}/${c}.md`;
    if (!existsSync(join(root, f))) throw new Error(`missing ${f}`);
    chapterFiles.push(f);
  }
  yamlPub.set(path, text);
}

const wt = mkdtempSync(join(tmpdir(), 'edumath-stage-'));
let ok = false;
try {
  git('worktree', 'add', '--detach', wt, 'HEAD');
  for (const f of chapterFiles) cpSync(join(root, f), join(wt, f));
  for (const [p, text] of yamlPub) writeFileSync(join(wt, p), text);
  symlinkSync(join(root, 'node_modules'), join(wt, 'node_modules'));
  if (existsSync(join(root, '.cache', 'tikz'))) cpSync(join(root, '.cache', 'tikz'), join(wt, '.cache', 'tikz'), { recursive: true });
  for (const cmd of [['run', 'check'], ['run', 'build']]) {
    try {
      execFileSync('npm', cmd, { cwd: wt, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    } catch (e) {
      const out = `${e.stdout ?? ''}\n${e.stderr ?? ''}`;
      const lines = out.split('\n').filter((l) => /error|Error|ERROR|unresolved|unknown|invalid|✗/.test(l));
      console.error(`npm ${cmd.join(' ')} failed:\n${lines.slice(0, 60).join('\n')}`);
      throw new Error('verification failed');
    }
  }
  ok = true;
} finally {
  try {
    git('worktree', 'remove', '--force', wt);
  } catch {
    rmSync(wt, { recursive: true, force: true });
    git('worktree', 'prune');
  }
}

if (ok) {
  git('add', ...chapterFiles);
  for (const [p, text] of yamlPub) {
    const blob = execFileSync('git', ['hash-object', '-w', '--stdin'], { cwd: root, input: text, encoding: 'utf8' }).trim();
    git('update-index', '--cacheinfo', `100644,${blob},${p}`);
    let work = readFileSync(join(root, p), 'utf8');
    for (const c of bySubject.get(p.split('/').pop().replace('.yaml', '')) ?? []) work = setStatus(work, c, 'reviewed');
    writeFileSync(join(root, p), work);
  }
  console.log(`staged ${chapterFiles.length} chapter(s): ${keys.join(', ')}`);
}
