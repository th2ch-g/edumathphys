import type { Root } from 'hast';
import { visit } from 'unist-util-visit';

/** Keep TeX as text and remove code wrappers before Expressive Code runs. */
export default function rehypeMathSource() {
  return (tree: Root): void => {
    visit(tree, 'element', (node) => {
      const code = node.tagName === 'pre' ? node.children[0] : node;
      if (code?.type !== 'element' || code.tagName !== 'code') return;

      const classes = code.properties.className;
      if (!Array.isArray(classes) || !classes.includes('language-math')) return;

      const display = node.tagName === 'pre' || classes.includes('math-display');
      node.tagName = display ? 'div' : 'span';
      node.properties = {
        className: [display ? 'math-display' : 'math-inline'],
        'data-pagefind-ignore': '',
      };
      node.children = code.children;
    });
  };
}
