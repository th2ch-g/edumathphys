import ELK from 'elkjs/lib/elk.bundled.js';

export interface GraphNode {
  id: string;
  label: string;
  url: string;
  className?: string;
}
export interface GraphEdge {
  from: string;
  to: string;
}
export interface LaidOutGraph {
  width: number;
  height: number;
  nodes: (GraphNode & { x: number; y: number; w: number; h: number })[];
  edges: { from: string; to: string; d: string }[];
}

const CHAR_W = 13;
const NODE_H = 34;

function textWidth(s: string): number {
  let w = 0;
  for (const ch of s) w += /[\x20-\x7e]/.test(ch) ? CHAR_W * 0.6 : CHAR_W;
  return Math.ceil(w);
}

/** Layered (Sugiyama-style) layout computed at build time. */
export async function layoutGraph(nodes: GraphNode[], edges: GraphEdge[], direction: 'RIGHT' | 'DOWN' = 'RIGHT'): Promise<LaidOutGraph> {
  const elk = new ELK();
  const ids = new Set(nodes.map((n) => n.id));
  const valid = edges.filter((e) => ids.has(e.from) && ids.has(e.to));
  const res = await elk.layout({
    id: 'root',
    layoutOptions: {
      'elk.algorithm': 'layered',
      'elk.direction': direction,
      'elk.layered.spacing.nodeNodeBetweenLayers': '56',
      'elk.spacing.nodeNode': '14',
      'elk.edgeRouting': 'POLYLINE',
      'elk.layered.nodePlacement.strategy': 'BRANDES_KOEPF',
    },
    children: nodes.map((n) => ({ id: n.id, width: textWidth(n.label) + 24, height: NODE_H })),
    edges: valid.map((e, i) => ({ id: `e${i}`, sources: [e.from], targets: [e.to] })),
  });
  const pos = new Map((res.children ?? []).map((c) => [c.id, c]));
  return {
    width: Math.ceil(res.width ?? 0),
    height: Math.ceil(res.height ?? 0),
    nodes: nodes.map((n) => {
      const p = pos.get(n.id)!;
      return { ...n, x: p.x ?? 0, y: p.y ?? 0, w: p.width ?? 0, h: p.height ?? 0 };
    }),
    edges: (res.edges ?? []).map((e, i) => {
      const sec = e.sections?.[0];
      const pts = sec ? [sec.startPoint, ...(sec.bendPoints ?? []), sec.endPoint] : [];
      return { from: valid[i].from, to: valid[i].to, d: smoothPath(pts) };
    }),
  };
}

function smoothPath(pts: { x: number; y: number }[]): string {
  if (pts.length === 0) return '';
  if (pts.length === 2) return `M${pts[0].x},${pts[0].y} L${pts[1].x},${pts[1].y}`;
  let d = `M${pts[0].x},${pts[0].y}`;
  for (let i = 1; i < pts.length - 1; i++) {
    const mid = { x: (pts[i].x + pts[i + 1].x) / 2, y: (pts[i].y + pts[i + 1].y) / 2 };
    d += ` Q${pts[i].x},${pts[i].y} ${i === pts.length - 2 ? `${pts[i + 1].x},${pts[i + 1].y}` : `${mid.x},${mid.y}`}`;
  }
  return d;
}
