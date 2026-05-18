// =============================================================
//  VariantScope — app.js
//  Deux backends selon l'assemblage :
//    • UCSC  : liftover + séquence  (canFam3/4/6 chien)
//    • Ensembl : liftover + séquence + variants  (tous les autres)
// =============================================================

// ----------- Configuration -----------------------------------

const SPECIES_CONFIG = {
    dog: {
        label:          'Chien (Canis lupus familiaris)',
        ensemblSpecies: 'canis_lupus_familiaris',
        genomes: [
            // Ensembl → pipeline complet (séquence + variants)
            { id: 'ROS_Cfam_1.0', label: 'ROS_Cfam_1.0 / Labrador SID07034 2020 (Ensembl — variants disponibles)',
              backend: 'ensembl', ensemblAsm: 'ROS_Cfam_1.0' },
            // UCSC → liftover vers GCF_014441545.1 (= ROS_Cfam_1.0) → variants Ensembl disponibles
            { id: 'canFam6', label: 'canFam6 / Dog10K_Boxer_Tasha',
              backend: 'ucsc', ucscDb: 'canFam6', ucscTarget: 'GCF_014441545.1' },
            { id: 'canFam5', label: 'canFam5 / UMICH_Zoey_3.1',
              backend: 'ucsc', ucscDb: 'canFam5', ucscTarget: 'GCF_014441545.1' },
            { id: 'canFam4', label: 'canFam4 / UU_Cfam_GSD_1.0',
              backend: 'ucsc', ucscDb: 'canFam4', ucscTarget: 'GCF_014441545.1' },
            { id: 'canFam3', label: 'canFam3 / CanFam3.1',
              backend: 'ucsc', ucscDb: 'canFam3', ucscTarget: 'GCF_014441545.1' }
        ],
        targetEnsemblAsm: 'ROS_Cfam_1.0'
    },
    cat: {
        label:          'Chat (Felis catus)',
        ensemblSpecies: 'felis_catus',
        genomes: [
            { id: 'Felis_catus_9.0', label: 'Felis_catus_9.0 / felCat9 (récent)',
              backend: 'ensembl', ensemblAsm: 'Felis_catus_9.0' },
            { id: 'Felis_catus_8.0', label: 'Felis_catus_8.0 / felCat8',
              backend: 'ensembl', ensemblAsm: 'Felis_catus_8.0' }
        ],
        targetEnsemblAsm: 'Felis_catus_9.0'
    },
    horse: {
        label:          'Cheval (Equus caballus)',
        ensemblSpecies: 'equus_caballus',
        genomes: [
            { id: 'EquCab3.0', label: 'EquCab3.0 / equCab3 (récent)',
              backend: 'ensembl', ensemblAsm: 'EquCab3.0' },
            { id: 'EquCab2.0', label: 'EquCab2.0 / equCab2',
              backend: 'ensembl', ensemblAsm: 'EquCab2.0' }
        ],
        targetEnsemblAsm: 'EquCab3.0'
    }
};

let _rawSeq    = '';
let _maskedSeq = '';

// ----------- Parsing de position ---------------------------------

function parsePosition(str) {
    str = str.trim();
    if (str.includes('_')) {
        const parts = str.split('_').map(Number);
        return { start: parts[0], end: parts[1] };
    }
    const n = Number(str);
    return { start: n, end: n };
}

// ----------- Formulaire --------------------------------------

function initForm() {
    const sel = document.getElementById('species');
    Object.entries(SPECIES_CONFIG).forEach(([key, cfg]) => {
        const opt = document.createElement('option');
        opt.value = key;
        opt.textContent = cfg.label;
        sel.appendChild(opt);
    });
    refreshGenomeOptions();
}

function refreshGenomeOptions() {
    const key = document.getElementById('species').value;
    const sel = document.getElementById('genome');
    sel.innerHTML = '';
    SPECIES_CONFIG[key].genomes.forEach(g => {
        const opt = document.createElement('option');
        opt.value = g.id;
        opt.textContent = g.label;
        sel.appendChild(opt);
    });
}

document.getElementById('species').addEventListener('change', refreshGenomeOptions);

// ----------- Soumission --------------------------------------

document.getElementById('variantForm').addEventListener('submit', async (e) => {
    e.preventDefault();

    const speciesKey = document.getElementById('species').value;
    const genomeId   = document.getElementById('genome').value;
    let   chrom      = document.getElementById('chromosome').value.trim();
    const posRange   = parsePosition(document.getElementById('position').value);
    const windowSize = Math.max(50, parseInt(document.getElementById('windowSize').value, 10) || 500);

    const species  = SPECIES_CONFIG[speciesKey];
    const genomeCfg = species.genomes.find(g => g.id === genomeId);

    const btn = document.getElementById('submitBtn');
    btn.disabled = true;
    resetUI();
    showEl('results');
    showEl('loading');
    hideEl('errorBox');

    try {
        if (genomeCfg.backend === 'ucsc') {
            await runUCSCPipeline(species, genomeCfg, chrom, posRange.start, posRange.end, windowSize);
        } else {
            await runEnsemblPipeline(species, genomeCfg, chrom, posRange.start, posRange.end, windowSize);
        }
    } catch (err) {
        showError(err.message || 'Erreur inattendue.');
    } finally {
        hideEl('loading');
        btn.disabled = false;
    }
});

// =============================================================
//  Pipeline UCSC (canFam3 / canFam4 / canFam6)
//  Séquence uniquement — pas de variants (coordonnées incompatibles Ensembl)
// =============================================================

async function runUCSCPipeline(species, genomeCfg, chrom, posStart, posEnd, windowSize) {
    // Normaliser "chr" (UCSC en a besoin)
    if (!/^chr/i.test(chrom)) chrom = 'chr' + chrom;

    let finalDb       = genomeCfg.ucscDb;
    let finalChrom    = chrom;
    let finalPosStart = posStart;
    let finalPosEnd   = posEnd;

    // LiftOver UCSC si l'assemblage n'est pas déjà la cible
    if (genomeCfg.ucscDb !== genomeCfg.ucscTarget) {
        setLoadingMsg('LiftOver UCSC en cours…');
        const lifted  = await ucscLiftOver(genomeCfg.ucscDb, genomeCfg.ucscTarget, chrom, posStart);
        finalDb       = genomeCfg.ucscTarget;
        finalChrom    = lifted.chrom;
        finalPosStart = lifted.position;
        finalPosEnd   = lifted.position + (posEnd - posStart);
        showLiftoverBanner(genomeCfg.ucscDb, genomeCfg.ucscTarget,
                           chrom, posStart, posEnd, finalChrom, finalPosStart, finalPosEnd);
    }

    // Fenêtre centrée sur le milieu de l'intervalle (UCSC : 0-based half-open)
    const midPos    = Math.round((finalPosStart + finalPosEnd) / 2);
    const winStart0 = Math.max(0, midPos - 1 - windowSize);
    const winEnd0   = midPos + windowSize;
    const winStart1 = winStart0 + 1;
    const mutIdxStart = finalPosStart - 1 - winStart0;
    const mutIdxEnd   = finalPosEnd   - 1 - winStart0;

    setLoadingMsg('Récupération de la séquence UCSC…');
    const sequence = await ucscSequence(finalDb, finalChrom, winStart0, winEnd0);

    _rawSeq = sequence;
    displaySequence(sequence, mutIdxStart, mutIdxEnd, finalChrom, winStart1, winEnd0);

    // GCF_014441545.1 = ROS_Cfam_1.0 → variants Ensembl sur les mêmes coordonnées
    if (finalDb === 'GCF_014441545.1') {
        setLoadingMsg('Récupération des variants Ensembl…');
        const chromEns = finalChrom.replace(/^chr/i, '');
        const variants = await ensemblVariants(species.ensemblSpecies, chromEns, winStart1, winEnd0);
        displayVariants(variants);
        const masked = buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
        _maskedSeq = masked.join('');
        displayMaskedSequence(masked, variants, mutIdxStart, mutIdxEnd, finalChrom, winStart1, winEnd0);
    } else {
        _maskedSeq = '';
        showVariantUnavailable(genomeCfg.id, 'ROS_Cfam_1.0');
    }
}

// =============================================================
//  Pipeline Ensembl (chat, cheval, chien en ROS_Cfam_1.0)
//  Séquence + variants complets
// =============================================================

async function runEnsemblPipeline(species, genomeCfg, chrom, posStart, posEnd, windowSize) {
    // Ensembl n'utilise pas le préfixe "chr"
    let chromEns = chrom.replace(/^chr/i, '');

    let finalAsm      = genomeCfg.ensemblAsm;
    let finalChrom    = chromEns;
    let finalPosStart = posStart;
    let finalPosEnd   = posEnd;

    // LiftOver Ensembl si nécessaire
    if (genomeCfg.ensemblAsm !== species.targetEnsemblAsm) {
        setLoadingMsg('LiftOver Ensembl en cours…');
        const lifted = await ensemblLiftOver(
            species.ensemblSpecies, genomeCfg.ensemblAsm,
            species.targetEnsemblAsm, chromEns, posStart
        );
        finalAsm      = species.targetEnsemblAsm;
        finalChrom    = lifted.chrom;
        finalPosStart = lifted.position;
        finalPosEnd   = lifted.position + (posEnd - posStart);
        showLiftoverBanner(genomeCfg.ensemblAsm, species.targetEnsemblAsm,
                           chromEns, posStart, posEnd, finalChrom, finalPosStart, finalPosEnd);
    }

    // Fenêtre centrée sur le milieu de l'intervalle (Ensembl : 1-based inclusif)
    const midPos    = Math.round((finalPosStart + finalPosEnd) / 2);
    const winStart1 = Math.max(1, midPos - windowSize);
    const winEnd1   = midPos + windowSize;
    const mutIdxStart = finalPosStart - winStart1;
    const mutIdxEnd   = finalPosEnd   - winStart1;

    setLoadingMsg('Récupération de la séquence Ensembl…');
    const sequence = await ensemblSequence(species.ensemblSpecies, finalChrom, winStart1, winEnd1);

    setLoadingMsg('Récupération des variants Ensembl…');
    const variants = await ensemblVariants(species.ensemblSpecies, finalChrom, winStart1, winEnd1);

    _rawSeq = sequence;
    displaySequence(sequence, mutIdxStart, mutIdxEnd, finalChrom, winStart1, winEnd1);
    displayVariants(variants);

    const masked = buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
    _maskedSeq = masked.join('');
    displayMaskedSequence(masked, variants, mutIdxStart, mutIdxEnd, finalChrom, winStart1, winEnd1);
}

// =============================================================
//  Fonctions API UCSC
// =============================================================

async function ucscLiftOver(fromDb, toDb, chrom, position) {
    // UCSC : coordonnées 0-based half-open pour une seule base
    // Appel via proxy Netlify pour éviter les restrictions CORS
    const start = position - 1;
    const end   = position;
    const url = `/.netlify/functions/liftover` +
                `?fromDb=${fromDb}&toDb=${toDb}` +
                `&chrom=${encodeURIComponent(chrom)}&start=${start}&end=${end}`;

    const resp = await fetch(url);
    if (!resp.ok) {
        let detail = '';
        try { const d = await resp.json(); detail = d.error || ''; } catch {}
        throw new Error(`LiftOver UCSC ${fromDb}→${toDb} : HTTP ${resp.status}${detail ? ' — ' + detail : ''}`);
    }
    const data = await resp.json();
    if (data.error) throw new Error('LiftOver UCSC : ' + data.error);
    if (!data.mappedCoordinates || data.mappedCoordinates.length === 0) {
        throw new Error(`LiftOver UCSC : aucune correspondance pour ${chrom}:${position} (${fromDb}→${toDb}).`);
    }

    const m = data.mappedCoordinates[0];
    return {
        chrom:    m.chrom,
        position: m.start + 1   // 0-based → 1-based
    };
}

async function ucscSequence(db, chrom, start0, end0) {
    // UCSC : 0-based half-open
    const url = `https://api.genome.ucsc.edu/getData/sequence` +
                `?genome=${db}&chrom=${encodeURIComponent(chrom)}&start=${start0}&end=${end0}`;

    const resp = await fetch(url);
    if (!resp.ok) {
        let detail = '';
        try { const d = await resp.json(); detail = d.error || ''; } catch {}
        throw new Error(`Séquence UCSC : HTTP ${resp.status}${detail ? ' — ' + detail : ''}`);
    }
    const data = await resp.json();
    if (data.error) throw new Error('Séquence UCSC : ' + data.error);
    return data.dna.toUpperCase();
}

// =============================================================
//  Fonctions API Ensembl
// =============================================================

const ENS_HEADERS = { 'Accept': 'application/json' };

async function ensemblFetch(url) {
    const resp = await fetch(url, { headers: ENS_HEADERS });
    if (!resp.ok) {
        let detail = '';
        try { const e = await resp.json(); detail = e.error || e.message || ''; } catch {}
        throw new Error(`Ensembl HTTP ${resp.status}${detail ? ' — ' + detail : ''}`);
    }
    return resp.json();
}

async function ensemblLiftOver(ensemblSpecies, fromAsm, toAsm, chrom, position) {
    const url = `https://rest.ensembl.org/map/${ensemblSpecies}` +
                `/${fromAsm}/${chrom}:${position}..${position}/${toAsm}`;
    const data = await ensemblFetch(url);

    if (!data.mappings || data.mappings.length === 0) {
        throw new Error(`LiftOver ${fromAsm}→${toAsm} : aucune correspondance pour ${chrom}:${position}.`);
    }
    const m = data.mappings[0].mapped;
    return { chrom: m.seq_region_name, position: m.start };
}

async function ensemblSequence(ensemblSpecies, chrom, start1, end1) {
    const url = `https://rest.ensembl.org/sequence/region/${ensemblSpecies}` +
                `/${chrom}:${start1}..${end1}?type=genomic`;
    const data = await ensemblFetch(url);
    return data.seq.toUpperCase();
}

async function ensemblVariants(ensemblSpecies, chrom, start1, end1) {
    const url = `https://rest.ensembl.org/overlap/region/${ensemblSpecies}` +
                `/${chrom}:${start1}-${end1}?feature=variation`;
    const resp = await fetch(url, { headers: ENS_HEADERS });
    if (resp.status === 404) return [];
    if (!resp.ok) {
        let detail = '';
        try { const e = await resp.json(); detail = e.error || ''; } catch {}
        throw new Error(`Variants Ensembl HTTP ${resp.status}${detail ? ' — ' + detail : ''}`);
    }
    const data = await resp.json();
    return Array.isArray(data) ? data : [];
}

// =============================================================
//  Séquence masquée
// =============================================================

function buildMaskedSeq(seq, variants, regionStart1, mutIdxStart, mutIdxEnd) {
    // Retourne un tableau de tokens :
    //   - position(s) ciblée(s) : [ref/var] (position unique) ou bases surlignées (intervalle)
    //   - autres variants       : N
    //   - reste                 : base originale
    const tokens = seq.split('');

    // 1. Pour une position unique : chercher le label allélique dans Ensembl
    let mutLabel = null;
    if (mutIdxStart === mutIdxEnd) {
        for (const v of variants) {
            const i0 = v.start - regionStart1;
            const i1 = (v.end !== undefined ? v.end : v.start) - regionStart1;
            if (mutIdxStart >= i0 && mutIdxStart <= i1) {
                const alleles = Array.isArray(v.alleles) ? v.alleles :
                                (typeof v.alleles === 'string' ? v.alleles.split('/') : []);
                if (alleles.length > 0) {
                    mutLabel = '[' + alleles.join('/') + ']';
                    break;
                }
            }
        }
    }

    // 2. Masquer les variants hors de l'intervalle ciblé
    variants.forEach(v => {
        const i0 = v.start - regionStart1;
        const i1 = (v.end !== undefined ? v.end : v.start) - regionStart1;
        for (let i = i0; i <= i1; i++) {
            if (i < 0 || i >= tokens.length) continue;
            if (i >= mutIdxStart && i <= mutIdxEnd) continue;
            tokens[i] = 'N';
        }
    });

    // 3. Annoter la/les position(s) ciblée(s)
    if (mutIdxStart === mutIdxEnd) {
        tokens[mutIdxStart] = mutLabel || '[' + seq[mutIdxStart] + '/?]';
    } else {
        // Intervalle : chercher un variant qui chevauche la plage
        let rangeLabel = null;
        for (const v of variants) {
            const i0 = v.start - regionStart1;
            const i1 = (v.end !== undefined ? v.end : v.start) - regionStart1;
            if (i0 <= mutIdxEnd && i1 >= mutIdxStart) {
                const alleles = Array.isArray(v.alleles) ? v.alleles :
                                (typeof v.alleles === 'string' ? v.alleles.split('/') : []);
                if (alleles.length > 0) { rangeLabel = '[' + alleles.join('/') + ']'; break; }
            }
        }
        const refBases = seq.slice(mutIdxStart, mutIdxEnd + 1);
        tokens[mutIdxStart] = rangeLabel || '[' + refBases + '/?]';
        for (let i = mutIdxStart + 1; i <= mutIdxEnd; i++) tokens[i] = '';
    }

    return tokens;
}

// =============================================================
//  Affichage
// =============================================================

function displaySequence(seq, mutIdxStart, mutIdxEnd, chrom, start1, end1) {
    document.getElementById('coordsDisplay').textContent =
        `${chrom}:${fmt(start1)}–${fmt(end1)}  ·  ${seq.length} pb`;
    const hMap = new Map();
    for (let i = mutIdxStart; i <= mutIdxEnd; i++) hMap.set(i, 'pos');
    document.getElementById('sequenceDisplay').innerHTML = renderSeq(seq, hMap);
    showEl('sequenceCard');
}

function displayVariants(variants) {
    document.getElementById('variantCount').textContent =
        `${variants.length} polymorphisme(s) trouvé(s) dans cette région`;

    const tableEl = document.getElementById('variantTable');
    if (variants.length === 0) {
        tableEl.innerHTML = '<p class="no-data">Aucun polymorphisme référencé dans cette région.</p>';
    } else {
        const rows = variants.slice(0, 500).map(v => {
            const name    = v.id || v.variation_name || '–';
            const link    = /^rs\d+/.test(name)
                            ? `<a href="https://www.ensembl.org/id/${name}" target="_blank" rel="noopener">${name}</a>`
                            : name;
            const alleles = Array.isArray(v.alleles) ? v.alleles.join('/') :
                            (typeof v.alleles === 'string' ? v.alleles : '–');
            const conseq  = Array.isArray(v.consequence_type) ? v.consequence_type[0]
                            : (v.consequence_type || v.feature_type || '–');
            return `<tr>
                <td>${link}</td>
                <td>${v.seq_region_name || '–'}</td>
                <td>${fmt(v.start) ?? '–'}</td>
                <td><span class="type-badge">${conseq}</span></td>
                <td>${alleles}</td>
            </tr>`;
        }).join('');

        tableEl.innerHTML = `
            <table>
                <thead><tr>
                    <th>ID</th><th>Chr</th><th>Position</th><th>Type</th><th>Allèles</th>
                </tr></thead>
                <tbody>${rows}</tbody>
            </table>
            ${variants.length > 500
              ? `<p class="note">Affichage limité aux 500 premiers sur ${variants.length} variants.</p>`
              : ''}`;
    }
    showEl('variantsCard');
}

function showVariantUnavailable(fromAsm, targetAsm) {
    document.getElementById('variantCount').textContent = 'Variants non disponibles';
    document.getElementById('variantTable').innerHTML =
        `<p class="no-data">Les variants Ensembl sont indexés en <strong>${targetAsm}</strong>. ` +
        `Avec l'assemblage <strong>${fromAsm}</strong>, les positions seraient incorrectes — ` +
        `ils ne sont donc pas affichés. Pour les variants, sélectionnez <strong>${targetAsm}</strong>.</p>`;
    showEl('variantsCard');
}

function displayMaskedSequence(tokens, variants, mutIdxStart, mutIdxEnd, chrom, start1, end1) {
    document.getElementById('maskedCoordsDisplay').textContent =
        `${chrom}:${fmt(start1)}–${fmt(end1)}  ·  ${variants.length} polymorphisme(s) annoté(s)`;

    const hMap = new Map();
    for (let i = 0; i < tokens.length; i++) {
        if (tokens[i] === 'N') hMap.set(i, 'n');
    }
    for (let i = mutIdxStart; i <= mutIdxEnd; i++) hMap.set(i, 'pos');

    document.getElementById('maskedDisplay').innerHTML = renderTokenSeq(tokens, hMap);
    showEl('maskedCard');
}

function renderTokenSeq(tokens, highlightMap) {
    const LINE = 60, BLOCK = 10;
    let html = '<div class="seq-display">';
    for (let ls = 0; ls < tokens.length; ls += LINE) {
        const le = Math.min(ls + LINE, tokens.length);
        html += `<div class="seq-line"><span class="seq-coord">${ls + 1}</span><span class="seq-bases">`;
        let i = ls;
        while (i < le) {
            const posInLine = i - ls;
            if (posInLine > 0 && posInLine % BLOCK === 0) html += ' ';
            const type = highlightMap.get(i);
            if (type === 'pos') {
                let run = tokens[i++];
                while (i < le && highlightMap.get(i) === 'pos') {
                    if ((i - ls) % BLOCK === 0) run += ' ';
                    run += tokens[i++];
                }
                html += `<span class="highlight-pos">${run}</span>`;
            } else if (type === 'n') {
                html += `<span class="highlight-n">${tokens[i++]}</span>`;
            } else {
                html += tokens[i++];
            }
        }
        html += '</span></div>';
    }
    return html + '</div>';
}

// =============================================================
//  Rendu séquence
// =============================================================

function renderSeq(seq, highlightMap) {
    const LINE = 60, BLOCK = 10;
    let html = '<div class="seq-display">';
    for (let ls = 0; ls < seq.length; ls += LINE) {
        const le = Math.min(ls + LINE, seq.length);
        html += `<div class="seq-line"><span class="seq-coord">${ls + 1}</span><span class="seq-bases">`;
        let i = ls;
        while (i < le) {
            const posInLine = i - ls;
            if (posInLine > 0 && posInLine % BLOCK === 0) html += ' ';
            const type = highlightMap.get(i);
            if (type === 'pos') {
                let run = seq[i++];
                while (i < le && highlightMap.get(i) === 'pos') {
                    if ((i - ls) % BLOCK === 0) run += ' ';
                    run += seq[i++];
                }
                html += `<span class="highlight-pos">${run}</span>`;
            } else if (type === 'n') {
                html += `<span class="highlight-n">${seq[i++]}</span>`;
            } else {
                html += seq[i++];
            }
        }
        html += '</span></div>';
    }
    return html + '</div>';
}

// =============================================================
//  Copie
// =============================================================

function copySeq(which) {
    const seq = which === 'raw' ? _rawSeq : _maskedSeq;
    if (!seq) return;
    navigator.clipboard.writeText(seq).then(() => {
        document.querySelectorAll('.copy-btn').forEach(btn => {
            if (btn.getAttribute('onclick')?.includes(which)) {
                btn.textContent = 'Copié !';
                setTimeout(() => { btn.textContent = 'Copier'; }, 1500);
            }
        });
    }).catch(() => alert('Copie automatique bloquée. Sélectionnez manuellement.'));
}

// =============================================================
//  Helpers UI
// =============================================================

function showEl(id) { document.getElementById(id)?.classList.remove('hidden'); }
function hideEl(id) { document.getElementById(id)?.classList.add('hidden'); }

function resetUI() {
    ['sequenceCard', 'variantsCard', 'maskedCard', 'liftoverInfo', 'loading'].forEach(hideEl);
}

function showError(msg) {
    const box = document.getElementById('errorBox');
    box.textContent = '⚠ ' + msg;
    box.classList.remove('hidden');
    hideEl('loading');
}

function setLoadingMsg(msg) {
    const el = document.getElementById('loadingMsg');
    if (el) el.textContent = msg;
}

function showLiftoverBanner(fromAsm, toAsm, fromChrom, fromStart, fromEnd, toChrom, toStart, toEnd) {
    const fromCoord = fromStart === fromEnd
        ? `${fromChrom}:${fmt(fromStart)}`
        : `${fromChrom}:${fmt(fromStart)}–${fmt(fromEnd)}`;
    const toCoord = toStart === toEnd
        ? `${toChrom}:${fmt(toStart)}`
        : `${toChrom}:${fmt(toStart)}–${fmt(toEnd)}`;
    const el = document.getElementById('liftoverInfo');
    el.innerHTML = `<strong>LiftOver effectué :</strong> ${fromAsm} ${fromCoord} → ${toAsm} ${toCoord}`;
    showEl('liftoverInfo');
}

function fmt(n) {
    return typeof n === 'number' ? n.toLocaleString('fr-FR') : n;
}

// =============================================================
//  Démarrage
// =============================================================
initForm();
