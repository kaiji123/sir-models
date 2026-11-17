// Static endpoint: /examples/<id>.R for every library model, generated at build time.
import { LIBRARY } from '../../lib/library.js';
import { generate } from '../../lib/generator.js';

export function getStaticPaths() {
  return LIBRARY.map((m) => ({ params: { id: m.id }, props: { model: m } }));
}

export function GET({ props }) {
  return new Response(generate(props.model).code, { headers: { 'Content-Type': 'text/plain; charset=utf-8' } });
}
