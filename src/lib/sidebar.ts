import { AREAS, LEVELS, loadCurriculum, type Area, type Level } from './curriculum.ts';

export const GUIDE_LINKS = [
  { label: 'このサイトについて', link: '/guide/' },
  { label: '学習ロードマップ', link: '/map/' },
  { label: '学習進捗', link: '/progress/' },
  { label: '執筆状況', link: '/coverage/' },
  { label: '用語索引', link: '/glossary/' },
  { label: '記号表', link: '/notation/' },
  { label: '記法と規約', link: '/guide/conventions/' },
];

export function buildSidebar() {
  const cur = loadCurriculum();
  const areaGroups = (Object.keys(AREAS) as Area[]).map((area) => ({
    label: AREAS[area],
    items: (Object.keys(LEVELS) as Level[])
      .map((level) => ({
        label: LEVELS[level],
        items: cur.subjects
          .filter((s) => s.area === area && s.level === level)
          .map((s) => ({
            label: s.title,
            collapsed: true,
            items: [
              { label: '概要', link: s.url },
              ...s.chapters.map((c) => ({
                label: `${c.number}. ${c.title}`,
                link: c.url,
                attrs: c.exists ? {} : { class: 'stub-link' },
              })),
            ],
          })),
      }))
      .filter((g) => g.items.length > 0),
  }));
  return [{ label: 'はじめに', items: GUIDE_LINKS }, ...areaGroups];
}
