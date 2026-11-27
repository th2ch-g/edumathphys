import { loadStore, updateStore } from './store.ts';

function renderReadButton(btn: HTMLButtonElement, read: boolean): void {
  btn.setAttribute('aria-pressed', String(read));
  btn.textContent = read ? '読了済み ✓' : '読了にする';
}

export function initProgress(): void {
  const store = loadStore();
  for (const btn of document.querySelectorAll<HTMLButtonElement>('button.read-toggle')) {
    const doc = btn.dataset.doc!;
    renderReadButton(btn, Boolean(store.read[doc]));
    btn.addEventListener('click', () => {
      const s = updateStore((st) => {
        if (st.read[doc]) delete st.read[doc];
        else st.read[doc] = Date.now();
      });
      renderReadButton(btn, Boolean(s.read[doc]));
    });
  }

  for (const ex of document.querySelectorAll<HTMLElement>('.exercise[data-exid]')) {
    const id = ex.dataset.exid!;
    const head = ex.querySelector('.env-head');
    if (!head) continue;
    const label = document.createElement('label');
    label.className = 'exercise-done';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.checked = Boolean(store.ex[id]);
    box.addEventListener('change', () => {
      updateStore((st) => {
        if (box.checked) st.ex[id] = Date.now();
        else delete st.ex[id];
      });
    });
    label.append(box, ' 解けた');
    head.append(label);
  }
}
