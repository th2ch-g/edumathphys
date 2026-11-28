import { katexOptions } from '../lib/katex-macros.ts';

export async function renderMath(root: ParentNode): Promise<void> {
  const sources = root.querySelectorAll<HTMLElement>('.math-inline, .math-display');
  if (sources.length === 0) return;

  const katex = await import('katex').then((module) => module.default, () => undefined);
  const output = document.createElement('span');
  for (const source of sources) {
    if (katex) {
      try {
        katex.render(source.textContent ?? '', output, {
          ...katexOptions,
          macros: { ...katexOptions.macros },
          displayMode: source.classList.contains('math-display'),
          throwOnError: false,
        });
        source.replaceWith(...output.childNodes);
        continue;
      } catch (error) {
        console.error('KaTeX rendering failed', error);
      }
    }
    source.setAttribute('data-math-fallback', '');
  }
}
