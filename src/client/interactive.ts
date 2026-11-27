type Mount = (el: HTMLElement, params: Record<string, string>) => void | Promise<void>;

const modules = import.meta.glob<{ default: Mount }>('../interactive/*.ts');

async function mount(fig: HTMLElement): Promise<void> {
  const name = fig.dataset.name!;
  const loader = modules[`../interactive/${name}.ts`];
  const canvas = fig.querySelector<HTMLElement>('.interactive-canvas');
  if (!loader || !canvas) return;
  const params = JSON.parse(fig.dataset.params ?? '{}') as Record<string, string>;
  canvas.replaceChildren();
  try {
    const mod = await loader();
    await mod.default(canvas, params);
  } catch (e) {
    canvas.textContent = `図の読み込みに失敗しました: ${(e as Error).message}`;
  }
}

export function initInteractive(): void {
  const figs = document.querySelectorAll<HTMLElement>('figure.interactive[data-name]');
  if (figs.length === 0) return;
  const io = new IntersectionObserver(
    (entries) => {
      for (const e of entries) {
        if (!e.isIntersecting) continue;
        io.unobserve(e.target);
        void mount(e.target as HTMLElement);
      }
    },
    { rootMargin: '200px' },
  );
  figs.forEach((f) => io.observe(f));
}
