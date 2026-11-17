// Static endpoint: /examples/<id>.json (the diagram in the editor's file format).
import { LIBRARY } from '../../lib/library.js';
import { serialize } from '../../lib/generator.js';

export function getStaticPaths() {
  return LIBRARY.map((m) => ({ params: { id: m.id }, props: { model: m } }));
}

export function GET({ props }) {
  return new Response(serialize(props.model), { headers: { 'Content-Type': 'application/json' } });
}
