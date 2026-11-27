import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import type { AstroIntegration } from 'astro';
import { loadCurriculum, loadReferences } from './curriculum.ts';
import { buildDocIndex, invalidateDocIndex } from './doc-index.ts';

export function collectErrors(): string[] {
  const cur = loadCurriculum({ fresh: true });
  const { refs, errors: refErrors } = loadReferences({ fresh: true });
  const index = buildDocIndex({ fresh: true });
  const subjectRefErrors = cur.subjects.flatMap((s) =>
    s.references.filter((k) => !refs.has(k)).map((k) => `${s.sourceFile}: unknown reference key "${k}"`),
  );
  return [...cur.errors, ...refErrors, ...subjectRefErrors, ...index.errors];
}

/** Validates the curriculum and content index before building, and refreshes it in dev. */
export function edumathIntegration(opts: { watch: string[] }): AstroIntegration {
  return {
    name: 'edumath',
    hooks: {
      'astro:config:setup': ({ command, logger }) => {
        const errors = collectErrors();
        if (errors.length === 0) return;
        const msg = `${errors.length} content error(s):\n- ${errors.join('\n- ')}`;
        if (command === 'build') throw new Error(msg);
        logger.warn(msg);
      },
      'astro:build:done': async ({ dir }) => {
        // Raster icon for platforms that do not support SVG favicons.
        const { default: sharp } = await import('sharp');
        const svg = readFileSync(join(process.cwd(), 'public', 'favicon.svg'));
        await sharp(svg, { density: 300 }).resize(180, 180).png().toFile(join(fileURLToPath(dir), 'apple-touch-icon.png'));
      },
      'astro:server:setup': ({ server, logger }) => {
        const onChange = (file: string) => {
          if (!opts.watch.some((d) => file.startsWith(d))) return;
          invalidateDocIndex();
          const errors = collectErrors();
          if (errors.length > 0) logger.warn(`${errors.length} content error(s):\n- ${errors.join('\n- ')}`);
        };
        server.watcher.on('change', onChange);
        server.watcher.on('add', onChange);
        server.watcher.on('unlink', onChange);
      },
    },
  };
}
