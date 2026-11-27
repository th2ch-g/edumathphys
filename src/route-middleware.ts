import { defineRouteMiddleware, type StarlightRouteData } from '@astrojs/starlight/route-data';

type Entry = StarlightRouteData['sidebar'][number];

const hasCurrent = (e: Entry): boolean => (e.type === 'link' ? e.isCurrent : e.entries.some(hasCurrent));

function sectionOf(pathname: string, base: string): string {
  const rel = pathname.startsWith(base) ? pathname.slice(base.length) : pathname;
  const parts = rel.split('/').filter(Boolean);
  return parts[0] === 'math' || parts[0] === 'physics' ? parts.slice(0, 2).join('/') : (parts[0] ?? '');
}

/** Keeps the sidebar small: only the current subject is expanded, other subjects become plain links. */
export const onRequest = defineRouteMiddleware((context) => {
  const route = context.locals.starlightRoute;
  const base = import.meta.env.BASE_URL.replace(/\/?$/, '/');
  const onSubjectPage = /^(math|physics)\//.test(sectionOf(context.url.pathname, base));

  route.sidebar = route.sidebar.map((top) => {
    if (top.type !== 'group') return top;
    const current = hasCurrent(top);
    const entries = top.entries.map((level) => {
      if (level.type !== 'group') return level;
      return {
        ...level,
        entries: level.entries.map((subject): Entry => {
          if (subject.type !== 'group') return subject;
          if (hasCurrent(subject)) return { ...subject, collapsed: false };
          const overview = subject.entries[0];
          const href = overview?.type === 'link' ? overview.href : '#';
          return { type: 'link', label: subject.label, href, isCurrent: false, badge: undefined, attrs: {} };
        }),
      };
    });
    const collapsed = top.label === 'はじめに' ? onSubjectPage : !current;
    return { ...top, entries, collapsed };
  });

  const here = sectionOf(context.url.pathname, base);
  const { pagination } = route;
  if (pagination.prev && sectionOf(pagination.prev.href, base) !== here) pagination.prev = undefined;
  if (pagination.next && sectionOf(pagination.next.href, base) !== here) pagination.next = undefined;
});
