import { initChem, createContext } from '../src/chem.js';
import { simulate } from '../src/simulator.js';
await initChem();
const mols = process.argv.slice(2).length ? process.argv.slice(2) : [
  'CCO', 'CC(C)=O', 'OC(=O)c1ccccc1', 'CCOC(C)=O', 'O=Cc1ccccc1', 'CC#N', 'Nc1ccccc1', 'Cc1ccccc1',
  'O=C1CCCCC1', 'O=[N+]([O-])c1ccccc1', 'CCCCCC', 'C=CCO', 'C#CCCC', 'CC(N)=O', 'CNC(C)=O', 'C/C=C/C', 'C/C=C\C',
  'CC(=O)OC(C)=O', 'CC(=O)Cl', 'CCOCC', 'COc1ccccc1', 'Cc1ccc(C)cc1', 'CC(C)(C)O', 'c1ccncc1', 'CS(C)=O', 'O=C1CCC1', 'CCN(CC)CC', 'FC(F)(F)c1ccccc1', 'CCCl', 'CC(=O)[O-].[Na+]'
];
let fail = 0;
for (const s of mols) {
  const ctx = createContext(s);
  try {
    const r = simulate(ctx);
    const top = r.bands.filter(b => !b.overtone).sort((a,b)=>b.height-a.height).slice(0,4).map(b=>`${b.center}(${b.height})`).join(' ');
    console.log(`${s.padEnd(24)} ${ctx.formula.padEnd(10)} IHD=${ctx.ihd} | ${r.groups.map(g=>g.id+(g.variants.length?'['+g.variants.join(';')+']':'')).join(', ')}\n   top: ${top}`);
  } catch (e) { fail++; console.error('FALLO', s, e); }
  finally { ctx.dispose(); }
}
process.exit(fail ? 1 : 0);
