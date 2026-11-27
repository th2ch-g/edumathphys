import JXG from 'jsxgraph';

export { JXG };

export interface Palette {
  fg: string;
  muted: string;
  accent: string;
  accent2: string;
  ok: string;
}

export function palette(): Palette {
  const cs = getComputedStyle(document.documentElement);
  const v = (name: string, fallback: string) => cs.getPropertyValue(name).trim() || fallback;
  return {
    fg: v('--sl-color-white', '#222'),
    muted: v('--sl-color-gray-4', '#888'),
    accent: v('--edu-thm', '#3b82f6'),
    accent2: v('--edu-ex', '#f59e0b'),
    ok: v('--edu-ok', '#16a34a'),
  };
}

let counter = 0;

/** Creates a JSXGraph board inside `el` with theme-aware defaults. */
export function createBoard(el: HTMLElement, bbox: [number, number, number, number], opts: { height?: number; axis?: boolean } = {}) {
  const box = document.createElement('div');
  box.id = `jxg-board-${++counter}`;
  box.className = 'jxgbox';
  box.style.width = '100%';
  box.style.height = `${opts.height ?? 360}px`;
  box.style.border = 'none';
  el.append(box);
  const p = palette();
  (JXG.Options.text as unknown as Record<string, unknown>).strokeColor = p.fg;
  const attrs: Record<string, unknown> = {
    boundingbox: bbox,
    axis: opts.axis ?? true,
    showCopyright: false,
    showNavigation: false,
    keepAspectRatio: true,
    pan: { enabled: false },
    zoom: { wheel: false, needShift: true },
    defaultAxes: {
      x: { strokeColor: p.muted, ticks: { strokeColor: p.muted, label: { color: p.muted } } },
      y: { strokeColor: p.muted, ticks: { strokeColor: p.muted, label: { color: p.muted } } },
    },
  };
  const board = JXG.JSXGraph.initBoard(box.id, attrs as never);
  return { board, palette: p };
}

export function controls(el: HTMLElement): HTMLDivElement {
  const div = document.createElement('div');
  div.className = 'interactive-controls';
  el.append(div);
  return div;
}

export function slider(parent: HTMLElement, label: string, min: number, max: number, step: number, value: number, onInput: (v: number) => void): HTMLInputElement {
  const wrap = document.createElement('label');
  const span = document.createElement('span');
  const input = document.createElement('input');
  input.type = 'range';
  input.min = String(min);
  input.max = String(max);
  input.step = String(step);
  input.value = String(value);
  const out = document.createElement('output');
  const update = () => {
    out.textContent = Number(input.value).toFixed(2);
    onInput(Number(input.value));
  };
  input.addEventListener('input', update);
  span.textContent = label;
  wrap.append(span, input, out);
  parent.append(wrap);
  out.textContent = Number(value).toFixed(2);
  return input;
}
