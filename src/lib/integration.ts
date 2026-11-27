import type { AstroIntegration } from 'astro';
import { loadCurriculum, loadReferences } from './curriculum.ts';
import { buildDocIndex, invalidateDocIndex } from './doc-index.ts';

export function collectErrors(): string[] {
  const cur = loadCurriculum({ fresh: true });
  const { errors: refErrors } = loadReferences();
  const index = buildDocIndex({ fresh: true });
  return [...cur.errors, ...refErrors, ...index.errors];
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
