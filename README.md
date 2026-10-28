# IRAnalizador

**En línea:** https://ir.quimicaeduc.org

Aplicación web (Node.js) para el **análisis de espectros infrarrojos (IR) de compuestos orgánicos**.

- **Simulación** del espectro IR a partir de una estructura dibujada en el editor **JSME**.
- **Consulta de bases de datos externas**: identifica el compuesto en **PubChem** (nombre, CAS, SMILES, InChIKey) y descarga sus espectros IR experimentales del **NIST Chemistry WebBook** (JCAMP-DX).
- **Análisis de espectros propios** (JCAMP-DX o CSV/TXT de dos columnas).
- **Detección de picos y asignación** automática de cada banda, con **explicación detallada**: tipo de vibración, grupo funcional, intervalo esperado, desviación y el *porqué* químico de su posición (conjugación, tensión de anillo, puentes de hidrógeno, efectos inductivos, fase, resonancia de Fermi, etc.).
- Diagnóstico por regiones (X–H, triples enlaces, dobles enlaces, huella dactilar), por ausencia de bandas y por índice de deficiencia de hidrógeno.

## Cómo funciona

| Componente | Descripción |
|---|---|
| `src/chem.js` | RDKit (WebAssembly) para leer la estructura, calcular fórmula, masa, IHD, anillos y búsqueda de subestructuras (SMARTS). |
| `src/groups.js` | Base de conocimiento: ~35 reglas de grupos funcionales con sus bandas (posición, anchura, intensidad), correcciones por conjugación, tensión de anillo, patrón de sustitución aromática, estereoquímica de alquenos, fase gas/condensada, y su explicación. Incluye la tabla de correlación general. |
| `src/simulator.js` | Construye el espectro sumando perfiles pseudo-Voigt (absorbancia normalizada → %T). |
| `src/analysis.js` | Lectura JCAMP-DX/CSV, normalización a malla de 2 cm⁻¹, detección de picos (Savitzky–Golay + prominencia), asignación pico ↔ banda esperada, similitud y diagnóstico sin estructura. |
| `src/explain.js` | Genera las explicaciones en español. |
| `src/external.js` | PubChem PUG-REST y NIST WebBook (con caché y cola de peticiones para respetar sus límites). |
| `public/` | Interfaz: JSME + Plotly. |

## Ejecución local

Requiere Node.js ≥ 20.

```bash
npm install
npm start          # http://localhost:3000
npm test           # prueba rápida del simulador con ~30 moléculas
```

Variables de entorno opcionales:

| Variable | Valor por defecto | Uso |
|---|---|---|
| `PORT` | `3000` | Puerto HTTP |
| `NIST_MIN_INTERVAL_MS` | `1500` | Intervalo mínimo entre peticiones a NIST |

## API

| Método | Ruta | Cuerpo / parámetros |
|---|---|---|
| GET | `/api/health` | — |
| POST | `/api/simulate` | `{ structure: "SMILES o molfile", phase: "condensada" \| "gas" }` |
| GET | `/api/search` | `?q=ethyl acetate&type=auto\|name\|cas\|smiles` |
| POST | `/api/analyze` | `{ source: { type: "nist", id: "C141786", index: 2 }, structure?: "SMILES", phase? }` o `{ source: { x:[…], y:[…], yunits } }` |
| POST | `/api/upload` | multipart: `file`, `structure?`, `phase?`, `yMode?` (`auto`, `T%`, `T`, `A`) |
| POST | `/api/depict` | `{ structure, atoms: [índices] }` → SVG con átomos resaltados |

## Despliegue en Dokploy (VPS)

El repositorio incluye un `Dockerfile` listo para producción (Node 24, usuario sin privilegios, *healthcheck* en `/api/health`).

1. En Dokploy: **Create Project** → *IRAnalizador* → **Create Service → Application**.
2. **Provider**: GitHub (o *Git* público) → repositorio `PGBeermann/IRAnalizador`, rama `main`.
3. **Build Type**: `Dockerfile` (ruta `./Dockerfile`, contexto `.`).
4. **Environment**: `PORT=3000` (opcional).
5. **Domains**: añada el dominio (en producción: `ir.quimicaeduc.org`), **Container Port** `3000`, HTTPS con Let's Encrypt. El registro DNS `A` del subdominio debe apuntar a la IP del VPS.
6. **Deploy**. Active *Autodeploy* para que cada `git push` a `main` vuelva a desplegar.

## Fuentes de datos y licencias

- Espectros experimentales: [NIST Chemistry WebBook, SRD 69](https://webbook.nist.gov/chemistry/) — © U.S. Secretary of Commerce. Se descargan bajo demanda y no se redistribuyen.
- Identificación de compuestos: [PubChem](https://pubchem.ncbi.nlm.nih.gov/) (NCBI).
- [JSME Molecule Editor](https://jsme-editor.github.io/) (BSD), [RDKit](https://www.rdkit.org/) (BSD), [Plotly.js](https://plotly.com/javascript/) (MIT).
- [jcampconverter](https://github.com/cheminfo/jcampconverter) (CC BY-NC-SA 4.0 — uso no comercial).

El código propio de IRAnalizador se distribuye con licencia MIT. Herramienta docente: el espectro simulado es una aproximación empírica por frecuencias de grupo y no sustituye a un espectro experimental ni a un cálculo mecanocuántico.
