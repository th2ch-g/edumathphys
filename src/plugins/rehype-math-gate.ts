import type { VFile } from 'vfile';

/** rehype-katex only reports parse errors as messages; turn them into build failures. */
export default function rehypeMathGate() {
  return function transformer(_tree: unknown, file: VFile): void {
    const bad = file.messages.filter((m) => m.source === 'rehype-katex');
    if (bad.length > 0) {
      const details = bad.map((m) => `${m.place ? `${JSON.stringify(m.place)} ` : ''}${m.cause instanceof Error ? m.cause.message : m.reason}`).join('\n');
      file.fail(`KaTeX could not render ${bad.length} formula(s):\n${details}`);
    }
  };
}
