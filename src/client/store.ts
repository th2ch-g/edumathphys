const KEY = 'edumathphys:v1';

export interface ProgressStore {
  read: Record<string, number>;
  quiz: Record<string, { ok: boolean; t: number }>;
  ex: Record<string, number>;
}

function empty(): ProgressStore {
  return { read: {}, quiz: {}, ex: {} };
}

export function loadStore(): ProgressStore {
  try {
    const raw = localStorage.getItem(KEY);
    if (!raw) return empty();
    const parsed = JSON.parse(raw) as Partial<ProgressStore>;
    return { read: parsed.read ?? {}, quiz: parsed.quiz ?? {}, ex: parsed.ex ?? {} };
  } catch {
    return empty();
  }
}

export function updateStore(fn: (s: ProgressStore) => void): ProgressStore {
  const s = loadStore();
  fn(s);
  try {
    localStorage.setItem(KEY, JSON.stringify(s));
  } catch {
    // Storage may be unavailable (private mode); progress is then kept for this page view only.
  }
  document.dispatchEvent(new CustomEvent('edumath:progress', { detail: s }));
  return s;
}

export function clearStore(): void {
  try {
    localStorage.removeItem(KEY);
  } catch {
    // ignore
  }
  document.dispatchEvent(new CustomEvent('edumath:progress', { detail: empty() }));
}
