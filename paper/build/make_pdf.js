// SCRIPT: IQGT-PDF-V4
const fs = require('fs'), path = require('path'); const { chromium } = require('playwright');
const F = path.resolve(__dirname, 'node_modules/@fontsource')  // npm i @fontsource/newsreader@5.1.0 @fontsource/spectral@5.1.0 @fontsource/jetbrains-mono@5.1.0 @fontsource/archivo-black@5.1.0;
const face = (fam, file, w, st) => `@font-face{font-family:'${fam}';src:url('file://${F}/${file}') format('woff2');font-weight:${w};font-style:${st};}`;
const faces = [face('Newsreader','newsreader/files/newsreader-latin-400-normal.woff2',400,'normal'), face('Newsreader','newsreader/files/newsreader-latin-500-normal.woff2',500,'normal'),
 face('Newsreader','newsreader/files/newsreader-latin-600-normal.woff2',600,'normal'), face('Newsreader','newsreader/files/newsreader-latin-400-italic.woff2',400,'italic'),
 face('Newsreader','newsreader/files/newsreader-latin-500-italic.woff2',500,'italic'), face('Spectral','spectral/files/spectral-latin-500-normal.woff2',500,'normal'),
 face('Spectral','spectral/files/spectral-latin-600-normal.woff2',600,'normal'), face('JetBrains Mono','jetbrains-mono/files/jetbrains-mono-latin-400-normal.woff2',400,'normal'),
 face('JetBrains Mono','jetbrains-mono/files/jetbrains-mono-latin-500-normal.woff2',500,'normal'), face('Archivo Black','archivo-black/files/archivo-black-latin-400-normal.woff2',400,'normal')].join('\n');
const print = `@page{size:Letter;margin:0.55in 0.6in;background:#070a12;}
html,body{background:#070a12 !important;-webkit-print-color-adjust:exact;print-color-adjust:exact;}
main.paper{max-width:none !important;margin:0 !important;padding:0 !important;}
figure,.thm,.proof,mjx-container[display="true"]{break-inside:avoid;}
h2,h3{break-after:avoid;}`;
(async () => { const [src, out] = process.argv.slice(2);
 let html = fs.readFileSync(src, 'utf8').replace(/<link rel="stylesheet" href="https:\/\/fonts\.googleapis\.com[^>]*>/, '').replace('</style>', faces + '\n' + print + '\n</style>');
 const tmp = path.resolve(path.dirname(out), '_print.html'); fs.writeFileSync(tmp, html);
 const b = await chromium.launch(); const p = await b.newPage(); await p.goto('file://' + tmp, { waitUntil: 'load' }); await p.evaluate(() => document.fonts.ready);
 await p.emulateMedia({ media: 'print' });
 await p.pdf({ path: out, format: 'Letter', printBackground: true, preferCSSPageSize: true });
 await b.close(); fs.unlinkSync(tmp); console.log('wrote', out); })();
