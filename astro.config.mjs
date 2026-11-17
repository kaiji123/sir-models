// @ts-check
import { defineConfig } from 'astro/config';

// Fully static output: `npm run build` produces dist/, deployable to any static host
// (Netlify, Vercel, Cloudflare Pages, GitHub Pages, S3...). Code generation runs in the browser.
export default defineConfig({
  output: 'static',
});
