/**
 * Fast consistency check without a full site build:
 * curriculum schema/DAG, content <-> curriculum, labels, refs, citations, terms, KaTeX syntax.
 */
import { readFileSync } from 'node:fs';
import { relative } from 'node:path';
import katex from 'katex';
import { visit } from 'unist-util-visit';
import { collectErrors } from '../src/lib/integration.ts';
import { buildDocIndex } from '../src/lib/doc-index.ts';
import { katexOptions } from '../src/lib/katex-macros.ts';
import { extractEquationLabel, parseMarkdown } from '../src/lib/markdown.ts';

const errors = collectErrors();
const { docs, notation } = buildDocIndex();

for (const doc of docs.values()) {
  const rel = relative(process.cwd(), doc.file);
  const { tree } = parseMarkdown(readFileSync(doc.file, 'utf8'));
  visit(tree, (node) => {
    if (node.type !== 'math' && node.type !== 'inlineMath') return;
    const value = node.type === 'math' ? extractEquationLabel(node.value).value : node.value;
    try {
      katex.renderToString(value, { ...katexOptions, macros: { ...katexOptions.macros }, displayMode: node.type === 'math', throwOnError: true });
    } catch (e) {
      errors.push(`${rel}:${node.position?.start.line ?? '?'}: ${(e as Error).message}`);
    }
  });
}
for (const n of notation) {
  try {
    katex.renderToString(n.symbol, { ...katexOptions, macros: { ...katexOptions.macros }, throwOnError: true });
  } catch (e) {
    errors.push(`notation "${n.symbol}" (${n.url}): ${(e as Error).message}`);
  }
}

if (errors.length > 0) {
  console.error(`${errors.length} error(s):\n- ${errors.join('\n- ')}`);
  process.exit(1);
}
console.log(`OK: ${docs.size} documents checked`);
