import { initInteractive } from './interactive.ts';
import { renderMath } from './math.ts';
import { initProgress } from './progress.ts';
import { initQuizzes } from './quiz.ts';
import { initXrefPreview } from './xref-preview.ts';

void renderMath(document);
initProgress();
initQuizzes();
initXrefPreview();
initInteractive();
