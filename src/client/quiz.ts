import { loadStore, updateStore } from './store.ts';

function parseNumber(s: string): number {
  return Number(s.trim().replace(/[，,]/g, '').replace(/[−–]/g, '-').replace(/[０-９．]/g, (c) => String.fromCharCode(c.charCodeAt(0) - 0xfee0)));
}

function initQuiz(quiz: HTMLElement): void {
  const type = quiz.dataset.type as 'single' | 'multiple' | 'numeric';
  const qid = quiz.dataset.qid!;
  const options = [...quiz.querySelectorAll<HTMLLIElement>('.quiz-option')];
  const explanation = quiz.querySelector<HTMLElement>('.quiz-explanation');
  let input: HTMLInputElement | undefined;

  if (type === 'numeric') {
    input = document.createElement('input');
    input.type = 'text';
    input.inputMode = 'decimal';
    input.className = 'quiz-input';
    input.setAttribute('aria-label', '数値で解答');
    const wrap = document.createElement('p');
    wrap.className = 'quiz-input-row';
    wrap.append('答: ', input);
    (explanation ?? null) ? quiz.insertBefore(wrap, explanation) : quiz.append(wrap);
  } else {
    for (const li of options) {
      li.tabIndex = 0;
      li.setAttribute('role', type === 'single' ? 'radio' : 'checkbox');
      li.setAttribute('aria-checked', 'false');
      const toggle = () => {
        if (quiz.classList.contains('answered')) return;
        const on = li.getAttribute('aria-checked') !== 'true';
        if (type === 'single') options.forEach((o) => o.setAttribute('aria-checked', 'false'));
        li.setAttribute('aria-checked', String(on));
      };
      li.addEventListener('click', (e) => {
        if ((e.target as HTMLElement).closest('a')) return;
        toggle();
      });
      li.addEventListener('keydown', (e) => {
        if (e.key === ' ' || e.key === 'Enter') {
          e.preventDefault();
          toggle();
        }
      });
    }
  }

  const actions = document.createElement('div');
  actions.className = 'quiz-actions';
  const check = document.createElement('button');
  check.type = 'button';
  check.textContent = '答え合わせ';
  const reset = document.createElement('button');
  reset.type = 'button';
  reset.textContent = 'もう一度';
  reset.hidden = true;
  const result = document.createElement('span');
  result.className = 'quiz-result';
  result.setAttribute('aria-live', 'polite');
  actions.append(check, reset, result);
  explanation ? quiz.insertBefore(actions, explanation) : quiz.append(actions);

  const prev = loadStore().quiz[qid];
  if (prev) result.textContent = prev.ok ? '（前回: 正解）' : '（前回: 不正解）';

  check.addEventListener('click', () => {
    let ok: boolean;
    if (type === 'numeric') {
      const v = parseNumber(input!.value);
      if (input!.value.trim() === '' || Number.isNaN(v)) {
        result.textContent = '数値を入力してください';
        return;
      }
      const ans = Number(quiz.dataset.answer);
      const tol = Number(quiz.dataset.tol ?? '0');
      ok = Math.abs(v - ans) <= Math.max(tol, 1e-12 * Math.abs(ans));
      input!.readOnly = true;
    } else {
      const correct = new Set((quiz.dataset.correct ?? '').split(',').filter(Boolean).map(Number));
      const chosen = options.flatMap((o, i) => (o.getAttribute('aria-checked') === 'true' ? [i] : []));
      if (chosen.length === 0) {
        result.textContent = '選択肢を選んでください';
        return;
      }
      ok = chosen.length === correct.size && chosen.every((i) => correct.has(i));
      options.forEach((o, i) => {
        if (correct.has(i)) o.classList.add('is-correct');
        else if (o.getAttribute('aria-checked') === 'true') o.classList.add('is-wrong');
      });
    }
    quiz.classList.add('answered', ok ? 'ok' : 'ng');
    result.textContent = ok ? '正解' : '不正解';
    if (explanation) explanation.hidden = false;
    check.hidden = true;
    reset.hidden = false;
    updateStore((s) => {
      s.quiz[qid] = { ok, t: Date.now() };
    });
  });

  reset.addEventListener('click', () => {
    quiz.classList.remove('answered', 'ok', 'ng');
    options.forEach((o) => {
      o.classList.remove('is-correct', 'is-wrong');
      o.setAttribute('aria-checked', 'false');
    });
    if (input) {
      input.readOnly = false;
      input.value = '';
    }
    if (explanation) explanation.hidden = true;
    result.textContent = '';
    check.hidden = false;
    reset.hidden = true;
  });
}

export function initQuizzes(): void {
  document.querySelectorAll<HTMLElement>('.quiz[data-qid]').forEach(initQuiz);
}
