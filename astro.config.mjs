import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import { unified } from '@astrojs/markdown-remark';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import starlightLinksValidator from 'starlight-links-validator';
import { CONTENT_DIR, CURRICULUM_DIR, loadCurriculum } from './src/lib/curriculum.ts';
import { katexOptions } from './src/lib/katex-macros.ts';
import { buildSidebar } from './src/lib/sidebar.ts';
import { edumathIntegration } from './src/lib/integration.ts';
import remarkEduMath from './src/plugins/remark-edumath.ts';
import rehypeMathGate from './src/plugins/rehype-math-gate.ts';

const base = process.env.SITE_BASE ?? '/';

// Pages generated outside the docs collection cannot be checked by the links validator.
const customPages = new Set(['/', '/map/', '/glossary/', '/notation/', '/coverage/', '/progress/']);
for (const s of loadCurriculum().subjects) {
  customPages.add(s.url);
  for (const c of s.chapters) if (!c.exists) customPages.add(c.url);
}
const isCustomPage = ({ link }) => {
  const path = link.split('#')[0].replace(new RegExp(`^${base.replace(/\/$/, '')}`), '') || '/';
  return customPages.has(path.endsWith('/') ? path : `${path}/`);
};

export default defineConfig({
  site: process.env.SITE_URL,
  base,
  trailingSlash: 'always',
  markdown: {
    processor: unified({
      remarkPlugins: [remarkMath, [remarkEduMath, { include: [CONTENT_DIR] }]],
      rehypePlugins: [[rehypeKatex, katexOptions], rehypeMathGate],
      smartypants: false,
    }),
  },
  integrations: [
    edumathIntegration({ watch: [CONTENT_DIR, CURRICULUM_DIR] }),
    starlight({
      title: 'EduMathPhys',
      description: '大学数学・物理を学部から大学院まで日本語で体系的に独学するためのサイト',
      defaultLocale: 'root',
      locales: { root: { label: '日本語', lang: 'ja' } },
      sidebar: buildSidebar(),
      routeMiddleware: './src/route-middleware.ts',
      customCss: ['katex/dist/katex.min.css', './src/styles/custom.css'],
      components: {
        PageTitle: './src/components/PageTitle.astro',
        MarkdownContent: './src/components/MarkdownContent.astro',
      },
      lastUpdated: false,
      pagination: true,
      tableOfContents: { minHeadingLevel: 2, maxHeadingLevel: 3 },
      plugins: [starlightLinksValidator({ errorOnRelativeLinks: false, errorOnLocalLinks: true, exclude: isCustomPage })],
    }),
  ],
});
