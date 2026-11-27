import { controls, createBoard, slider } from './_board.ts';

/** Shows how a real 2x2 matrix maps the unit circle and a movable vector; eigen-directions are drawn when real. */
export default function mount(el: HTMLElement, params: Record<string, string>): void {
  const m = { a: Number(params.a ?? 2), b: Number(params.b ?? 1), c: Number(params.c ?? 1), d: Number(params.d ?? 2) };
  const { board, palette: p } = createBoard(el, [-4, 4, 4, -4]);
  const ctl = controls(el);

  board.create('curve', [(t: number) => Math.cos(t), (t: number) => Math.sin(t), 0, 2 * Math.PI], { strokeColor: p.muted, dash: 2, strokeWidth: 1 });
  board.create('curve', [(t: number) => m.a * Math.cos(t) + m.b * Math.sin(t), (t: number) => m.c * Math.cos(t) + m.d * Math.sin(t), 0, 2 * Math.PI], { strokeColor: p.accent2, strokeWidth: 2 });

  const v = board.create('glider', [Math.cos(0.6), Math.sin(0.6), board.create('circle', [[0, 0], 1], { visible: false })], { name: 'v', color: p.accent, size: 4 });
  const Av = board.create('point', [() => m.a * v.X() + m.b * v.Y(), () => m.c * v.X() + m.d * v.Y()], { name: 'Av', color: p.accent2, size: 4, fixed: true });
  board.create('arrow', [[0, 0], v], { strokeColor: p.accent, strokeWidth: 2 });
  board.create('arrow', [[0, 0], Av], { strokeColor: p.accent2, strokeWidth: 2 });

  const eig = () => {
    const tr = m.a + m.d;
    const det = m.a * m.d - m.b * m.c;
    const disc = tr * tr - 4 * det;
    if (disc < 0) return [];
    return [(tr + Math.sqrt(disc)) / 2, (tr - Math.sqrt(disc)) / 2].map((l) => {
      const [x, y] = Math.abs(m.b) > 1e-12 ? [m.b, l - m.a] : Math.abs(m.c) > 1e-12 ? [l - m.d, m.c] : l === m.a ? [1, 0] : [0, 1];
      const n = Math.hypot(x, y);
      return { l, x: x / n, y: y / n };
    });
  };
  for (const i of [0, 1]) {
    board.create('line', [[0, 0], [() => eig()[i]?.x ?? 0, () => eig()[i]?.y ?? 0]], {
      strokeColor: p.ok,
      strokeWidth: 1,
      dash: 1,
      visible: () => eig().length > i,
      fixed: true,
    });
  }
  const info = board.create('text', [-3.8, 3.6, () => {
    const e = eig();
    return e.length ? `固有値 λ₁ = ${e[0].l.toFixed(2)}, λ₂ = ${e[1].l.toFixed(2)}` : '実固有値なし（回転を含む）';
  }], { fixed: true, color: p.fg, fontSize: 13 });
  void info;

  for (const key of ['a', 'b', 'c', 'd'] as const) {
    slider(ctl, key, -3, 3, 0.1, m[key], (val) => {
      m[key] = val;
      board.update();
    });
  }
}
