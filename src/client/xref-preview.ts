import { renderMath } from './math.ts';

const pageCache = new Map<string, Promise<Document | undefined>>();

function fetchPage(url: string): Promise<Document | undefined> {
  let p = pageCache.get(url);
  if (!p) {
    p = fetch(url)
      .then((r) => (r.ok ? r.text() : Promise.reject(new Error(String(r.status)))))
      .then((html) => new DOMParser().parseFromString(html, 'text/html'))
      .catch(() => undefined);
    pageCache.set(url, p);
  }
  return p;
}

async function findTarget(link: HTMLAnchorElement): Promise<HTMLElement | undefined> {
  const url = new URL(link.href, location.href);
  const id = decodeURIComponent(url.hash.slice(1));
  if (!id) return undefined;
  if (url.pathname === location.pathname) return document.getElementById(id) ?? undefined;
  const doc = await fetchPage(url.pathname);
  return doc?.getElementById(id) ?? undefined;
}

export function initXrefPreview(): void {
  const pop = document.createElement('div');
  pop.className = 'xref-popover';
  pop.setAttribute('role', 'tooltip');
  pop.hidden = true;
  document.body.append(pop);

  let showTimer: number | undefined;
  let hideTimer: number | undefined;
  let active: HTMLAnchorElement | undefined;

  const hide = () => {
    pop.hidden = true;
    active?.removeAttribute('aria-describedby');
    active = undefined;
  };
  const scheduleHide = () => {
    window.clearTimeout(showTimer);
    hideTimer = window.setTimeout(hide, 200);
  };

  const show = async (link: HTMLAnchorElement) => {
    const target = await findTarget(link);
    if (!target || active !== link) return;
    const clone = target.cloneNode(true) as HTMLElement;
    clone.removeAttribute('id');
    clone.querySelectorAll('[id]').forEach((el) => el.removeAttribute('id'));
    clone.querySelectorAll('details').forEach((d) => d.remove());
    await renderMath(clone);
    if (active !== link) return;
    pop.replaceChildren(clone);
    pop.id = 'xref-popover';
    link.setAttribute('aria-describedby', pop.id);
    pop.hidden = false;
    const r = link.getBoundingClientRect();
    const width = Math.min(560, document.documentElement.clientWidth - 32);
    pop.style.width = `${width}px`;
    const left = Math.max(16, Math.min(r.left + window.scrollX, window.scrollX + document.documentElement.clientWidth - width - 16));
    const below = r.bottom + 8 + pop.offsetHeight < window.innerHeight;
    pop.style.left = `${left}px`;
    pop.style.top = below ? `${r.bottom + window.scrollY + 8}px` : `${r.top + window.scrollY - pop.offsetHeight - 8}px`;
  };

  const onEnter = (e: Event) => {
    const link = (e.target as HTMLElement).closest?.('a.xref') as HTMLAnchorElement | null;
    if (!link) return;
    window.clearTimeout(hideTimer);
    window.clearTimeout(showTimer);
    active = link;
    showTimer = window.setTimeout(() => void show(link), 250);
  };
  const onLeave = (e: Event) => {
    if ((e.target as HTMLElement).closest?.('a.xref')) scheduleHide();
  };

  document.addEventListener('mouseover', onEnter);
  document.addEventListener('focusin', onEnter);
  document.addEventListener('mouseout', onLeave);
  document.addEventListener('focusout', onLeave);
  pop.addEventListener('mouseenter', () => window.clearTimeout(hideTimer));
  pop.addEventListener('mouseleave', scheduleHide);
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') hide();
  });
}
