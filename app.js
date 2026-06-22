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
            // felCat9 = assemblage cible : séquence + variants EVA natifs
            { id: 'felCat9', label: 'felCat9 / Felis_catus_9.0',
              backend: 'ucsc', ucscDb: 'felCat9', variantTrack: 'evaSnp8' },
            { id: 'felCat8', label: 'felCat8 / Felis_catus_8.0 → liftover felCat9',
              backend: 'ucsc', ucscDb: 'felCat8', ucscTarget: 'felCat9', variantTrack: 'evaSnp8' },
            // Fca126 : liftover felCat9 non disponible → assemblage Ensembl natif
            { id: 'F.catus_Fca126_mat1.0',
              label: 'F.catus_Fca126_mat1.0 (Ensembl — variants Ensembl, pas de liftover felCat9)',
              backend: 'ensembl', ensemblAsm: 'F.catus_Fca126_mat1.0' }
        ],
        targetEnsemblAsm: 'F.catus_Fca126_mat1.0'
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

// Correspondance accessions RefSeq → noms karyotypiques pour les assemblages GCF
// (les fichiers chain UCSC utilisent les accessions NC_ comme noms de chromosomes)
const CHROM_ALIAS = {
    'GCF_014441545.1': {  // ROS_Cfam_1.0 (chien)
        'NC_051805.1': '1',  'NC_051806.1': '2',  'NC_051807.1': '3',
        'NC_051808.1': '4',  'NC_051809.1': '5',  'NC_051810.1': '6',
        'NC_051811.1': '7',  'NC_051812.1': '8',  'NC_051813.1': '9',
        'NC_051814.1': '10', 'NC_051815.1': '11', 'NC_051816.1': '12',
        'NC_051817.1': '13', 'NC_051818.1': '14', 'NC_051819.1': '15',
        'NC_051820.1': '16', 'NC_051821.1': '17', 'NC_051822.1': '18',
        'NC_051823.1': '19', 'NC_051824.1': '20', 'NC_051825.1': '21',
        'NC_051826.1': '22', 'NC_051827.1': '23', 'NC_051828.1': '24',
        'NC_051829.1': '25', 'NC_051830.1': '26', 'NC_051831.1': '27',
        'NC_051832.1': '28', 'NC_051833.1': '29', 'NC_051834.1': '30',
        'NC_051835.1': '31', 'NC_051836.1': '32', 'NC_051837.1': '33',
        'NC_051838.1': '34', 'NC_051839.1': '35', 'NC_051840.1': '36',
        'NC_051841.1': '37', 'NC_051842.1': '38', 'NC_051843.1': 'X'
    },
    'GCF_018350175.1': {  // F.catus_Fca126_mat1.0 (chat)
        'NC_001700.1': 'MT',
        'NC_058368.1': 'A1',  'NC_058369.1': 'A2',  'NC_058370.1': 'A3',
        'NC_058371.1': 'B1',  'NC_058372.1': 'B2',  'NC_058373.1': 'B3',
        'NC_058374.1': 'B4',  'NC_058375.1': 'C1',  'NC_058376.1': 'C2',
        'NC_058377.1': 'D1',  'NC_058378.1': 'D2',  'NC_058379.1': 'D3',
        'NC_058380.1': 'D4',  'NC_058381.1': 'E1',  'NC_058382.1': 'E2',
        'NC_058383.1': 'E3',  'NC_058384.1': 'F1',  'NC_058385.1': 'F2',
        'NC_058386.1': 'X'
    }
};

// Convertit un accession NC_ en nom de chromosome standard (sans préfixe)
function resolveChromName(db, chrom) {
    const map = CHROM_ALIAS[db];
    if (!map) return chrom;
    return map[chrom] || chrom;
}

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

    // LiftOver UCSC vers la cible (ex. canFam3 → GCF_014441545.1)
    if (genomeCfg.ucscTarget && genomeCfg.ucscDb !== genomeCfg.ucscTarget) {
        setLoadingMsg('LiftOver UCSC en cours…');
        const lifted  = await ucscLiftOver(genomeCfg.ucscDb, genomeCfg.ucscTarget, chrom, posStart);
        finalDb       = genomeCfg.ucscTarget;
        finalChrom    = lifted.chrom;
        finalPosStart = lifted.position;
        finalPosEnd   = lifted.position + (posEnd - posStart);
        showLiftoverBanner(genomeCfg.ucscDb, genomeCfg.ucscTarget,
                           chrom, posStart, posEnd, finalChrom, finalPosStart, finalPosEnd);
    }

    // Résoudre le nom de chromosome (accession NC_ → nom karyotypique si nécessaire)
    const baseChrom = resolveChromName(finalDb, finalChrom);
    const ucscChrom = /^NC_/.test(finalChrom) ? 'chr' + baseChrom : finalChrom;

    // Fenêtre centrée sur le milieu de l'intervalle (UCSC : 0-based half-open)
    const midPos    = Math.round((finalPosStart + finalPosEnd) / 2);
    const winStart0 = Math.max(0, midPos - 1 - windowSize);
    const winEnd0   = midPos + windowSize;
    const winStart1 = winStart0 + 1;
    const mutIdxStart = finalPosStart - 1 - winStart0;
    const mutIdxEnd   = finalPosEnd   - 1 - winStart0;

    setLoadingMsg('Récupération de la séquence UCSC…');
    const sequence = await ucscSequence(finalDb, ucscChrom, winStart0, winEnd0);

    _rawSeq = sequence;
    displaySequence(sequence, mutIdxStart, mutIdxEnd, ucscChrom, winStart1, winEnd0);

    // --- Cas 1 : cible GCF connue → variants Ensembl en coordonnées finales (chien)
    const GCF_ENSEMBL_MAP = { 'GCF_014441545.1': 'ROS_Cfam_1.0' };
    if (finalDb in GCF_ENSEMBL_MAP) {
        setLoadingMsg('Récupération des variants Ensembl…');
        const ensChrom = baseChrom.replace(/^chr/i, '');  // Ensembl n'accepte pas le préfixe "chr"
        const variants = await ensemblVariants(species.ensemblSpecies, ensChrom, winStart1, winEnd0);
        displayVariants(variants);
        const { tokens: masked, mutIdxStart: mStart, mutIdxEnd: mEnd } =
            buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
        _maskedSeq = masked.join('');
        displayMaskedSequence(masked, variants, mStart, mEnd, ucscChrom, winStart1, winEnd0);

    // --- Cas 2 : variants natifs UCSC (EVA SNP track) — toujours sur finalDb (felCat9)
    } else if (genomeCfg.variantTrack) {
        setLoadingMsg('Récupération des variants EVA…');
        const variants = await ucscVariants(
            finalDb, genomeCfg.variantTrack, ucscChrom, winStart0, winEnd0);
        displayVariants(variants);
        const { tokens: masked, mutIdxStart: mStart, mutIdxEnd: mEnd } =
            buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
        _maskedSeq = masked.join('');
        displayMaskedSequence(masked, variants, mStart, mEnd, ucscChrom, winStart1, winEnd0);

    } else {
        _maskedSeq = '';
        showVariantUnavailable(genomeCfg.id, species.targetEnsemblAsm);
    }

    // Position canFam3 (chien uniquement)
    if (species.ensemblSpecies === 'canis_lupus_familiaris') {
        if (genomeCfg.ucscDb === 'canFam3') {
            displayCanFam3Position(chrom, posStart, posEnd);
        } else {
            try {
                setLoadingMsg('Calcul de la position canFam3…');
                const cf3 = await ucscLiftOver(genomeCfg.ucscDb, 'canFam3', chrom, posStart);
                displayCanFam3Position(cf3.chrom, cf3.position, cf3.position + (posEnd - posStart));
            } catch { displayCanFam3Position(null, null, null); }
        }
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

    const { tokens: masked, mutIdxStart: mStart, mutIdxEnd: mEnd } =
        buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
    _maskedSeq = masked.join('');
    displayMaskedSequence(masked, variants, mStart, mEnd, finalChrom, winStart1, winEnd1);

    // Position canFam3 (chien uniquement)
    if (species.ensemblSpecies === 'canis_lupus_familiaris') {
        try {
            setLoadingMsg('Calcul de la position canFam3…');
            const cf3Chrom = 'chr' + finalChrom;
            const cf3 = await ucscLiftOver('GCF_014441545.1', 'canFam3', cf3Chrom, finalPosStart);
            displayCanFam3Position(cf3.chrom, cf3.position, cf3.position + (posEnd - posStart));
        } catch { displayCanFam3Position(null, null, null); }
    }
}

// =============================================================
//  Fonctions API UCSC
// =============================================================

async function ucscLiftOver(fromDb, toDb, chrom, position) {
    // Vérifier le cache localStorage avant d'appeler Netlify
    const cacheKey = `lo:${fromDb}:${toDb}:${chrom}:${position}`;
    try {
        const hit = localStorage.getItem(cacheKey);
        if (hit) return JSON.parse(hit);
    } catch {}

    const start = position - 1;
    const end   = position;
    const url = `/api/liftover` +
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
    const result = { chrom: m.chrom, position: m.start + 1 };

    try { localStorage.setItem(cacheKey, JSON.stringify(result)); } catch {}

    return result;
}

async function ucscVariants(db, track, chrom, start0, end0) {
    const url = `https://api.genome.ucsc.edu/getData/track` +
                `?genome=${db}&track=${track}` +
                `&chrom=${encodeURIComponent(chrom)}&start=${start0}&end=${end0}`;
    const resp = await fetch(url);
    if (!resp.ok) return [];
    const data = await resp.json();
    const items = data[track];
    if (!Array.isArray(items)) return [];

    // Convertir format UCSC EVA → format interne (compatible Ensembl)
    return items.map(v => ({
        id:               v.name,
        variation_name:   v.name,
        start:            v.chromStart + 1,   // 0-based → 1-based
        end:              v.chromEnd,          // 0-based end = 1-based end (insertion si end < start)
        alleles:          [v.ref || '-', v.alt || '-'],
        consequence_type: v.ucscClass || v.varClass || '–',
        seq_region_name:  v.chrom.replace(/^chr/i, '')
    }));
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
    // Retourne { tokens, mutIdxStart, mutIdxEnd } avec indices ajustés après insertions.
    // - position(s) ciblée(s) : [ref/var] ou bases surlignées
    // - insertions             : N inséré entre les deux bases flanquantes
    // - autres variants        : base remplacée par N
    const tokens = seq.split('');
    let adjMutStart = mutIdxStart;
    let adjMutEnd   = mutIdxEnd;

    // 1. Label allélique (indices originaux, avant toute modification)
    let mutLabel = null;
    if (mutIdxStart === mutIdxEnd) {
        for (const v of variants) {
            const i0 = v.start - regionStart1;
            const i1 = (v.end !== undefined ? v.end : v.start) - regionStart1;
            const covers = i1 < i0 ? (mutIdxStart === i0) : (mutIdxStart >= i0 && mutIdxStart <= i1);
            if (covers) {
                const alleles = Array.isArray(v.alleles) ? v.alleles :
                                (typeof v.alleles === 'string' ? v.alleles.split('/') : []);
                if (alleles.length > 0) { mutLabel = '[' + alleles.join('/') + ']'; break; }
            }
        }
    }

    // 2. Masquage — triés par position croissante pour suivre l'offset des insertions
    const sorted = [...variants].sort((a, b) => a.start - b.start);
    let offset = 0;

    for (const v of sorted) {
        const rawI0 = v.start - regionStart1;
        const rawI1 = (v.end !== undefined ? v.end : v.start) - regionStart1;

        if (rawI1 < rawI0) {
            // Insertion Ensembl (end < start) : insérer un N entre les bases flanquantes
            const insertPos = rawI0 + offset;
            if (insertPos < 0 || insertPos > tokens.length) continue;
            if (insertPos >= adjMutStart && insertPos <= adjMutEnd + 1) continue; // dans la plage ciblée
            tokens.splice(insertPos, 0, 'N');
            offset++;
            if (insertPos <= adjMutStart)      { adjMutStart++; adjMutEnd++; }
            else if (insertPos <= adjMutEnd)   { adjMutEnd++; }
        } else {
            const i0 = rawI0 + offset;
            const i1 = rawI1 + offset;
            for (let i = i0; i <= i1; i++) {
                if (i < 0 || i >= tokens.length) continue;
                if (i >= adjMutStart && i <= adjMutEnd) continue;
                tokens[i] = 'N';
            }
        }
    }

    // 3. Annoter la/les position(s) ciblée(s) (indices ajustés)
    if (mutIdxStart === mutIdxEnd) {
        tokens[adjMutStart] = mutLabel || '[' + seq[mutIdxStart] + '/?]';
    } else {
        let rangeLabel = null;
        for (const v of variants) {
            const i0 = v.start - regionStart1;
            const i1 = (v.end !== undefined ? v.end : v.start) - regionStart1;
            const effMin = Math.min(i0, i1);
            const effMax = Math.max(i0, i1);
            if (effMin <= mutIdxEnd && effMax >= mutIdxStart) {
                const alleles = Array.isArray(v.alleles) ? v.alleles :
                                (typeof v.alleles === 'string' ? v.alleles.split('/') : []);
                if (alleles.length > 0) { rangeLabel = '[' + alleles.join('/') + ']'; break; }
            }
        }
        const refBases = seq.slice(mutIdxStart, mutIdxEnd + 1);
        tokens[adjMutStart] = rangeLabel || '[' + refBases + '/?]';
        for (let i = adjMutStart + 1; i <= adjMutEnd; i++) tokens[i] = '';
    }

    return { tokens, mutIdxStart: adjMutStart, mutIdxEnd: adjMutEnd };
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

    const gcCount = (seq.match(/[GC]/g) || []).length;
    document.getElementById('gcContent').textContent =
        `GC : ${((gcCount / seq.length) * 100).toFixed(1)} %`;

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

function displayCanFam3Position(chrom, posStart, posEnd) {
    const el = document.getElementById('canfam3Content');
    if (!chrom) {
        el.textContent = 'Position canFam3 non disponible pour cette région.';
    } else {
        el.textContent = posStart === posEnd
            ? `${chrom}:${fmt(posStart)}`
            : `${chrom}:${fmt(posStart)}–${fmt(posEnd)}`;
    }
    showEl('canfam3Card');
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
    ['sequenceCard', 'variantsCard', 'maskedCard', 'liftoverInfo', 'loading', 'canfam3Card'].forEach(hideEl);
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
//  Traitement par lots
// =============================================================

// Correspondance noms de génome (colonne A du xlsx) → species/genomeId VariantScope
const GENOME_ALIAS = {
    'canfam3.1':             { species: 'dog',   genomeId: 'canFam3' },
    'canfam3':               { species: 'dog',   genomeId: 'canFam3' },
    'canfam4':               { species: 'dog',   genomeId: 'canFam4' },
    'canfam5':               { species: 'dog',   genomeId: 'canFam5' },
    'canfam6':               { species: 'dog',   genomeId: 'canFam6' },
    'ros_cfam_1.0':          { species: 'dog',   genomeId: 'ROS_Cfam_1.0' },
    'felcat9':               { species: 'cat',   genomeId: 'felCat9' },
    'felis_catus_9.0':       { species: 'cat',   genomeId: 'felCat9' },
    'felcat8':               { species: 'cat',   genomeId: 'felCat8' },
    'felis_catus_8.0':       { species: 'cat',   genomeId: 'felCat8' },
    'f.catus_fca126_mat1.0': { species: 'cat',   genomeId: 'F.catus_Fca126_mat1.0' },
    'equcab3.0':             { species: 'horse', genomeId: 'EquCab3.0' },
    'equcab3':               { species: 'horse', genomeId: 'EquCab3.0' },
    'equcab2.0':             { species: 'horse', genomeId: 'EquCab2.0' },
    'equcab2':               { species: 'horse', genomeId: 'EquCab2.0' },
};

let _batchOutputBuffer = null;

document.getElementById('batchFile').addEventListener('change', function () {
    const file = this.files[0];
    if (file) {
        document.getElementById('batchFileName').textContent = file.name;
        document.getElementById('batchRunBtn').disabled = false;
        _batchOutputBuffer = null;
        hideEl('batchDownload');
        hideEl('batchErrors');
        hideEl('batchProgress');
    }
});

async function runBatch() {
    const file = document.getElementById('batchFile').files[0];
    if (!file) return;

    document.getElementById('batchRunBtn').disabled = true;
    document.getElementById('batchErrors').innerHTML = '';
    hideEl('batchErrors');
    hideEl('batchDownload');
    showEl('batchProgress');
    const bar = document.getElementById('batchProgressBar');
    const msg = document.getElementById('batchProgressMsg');
    bar.style.width = '0%';
    bar.style.background = 'var(--primary)';
    msg.textContent = 'Chargement du fichier…';

    try {
        const arrayBuffer = await file.arrayBuffer();
        const workbook = new ExcelJS.Workbook();
        await workbook.xlsx.load(arrayBuffer);
        const ws = workbook.worksheets[0];

        // Collecter les lignes de données (col A = génome reconnu)
        const dataRows = [];
        ws.eachRow((row, rowNum) => {
            if (rowNum < 8) return;
            const genomeRaw = row.getCell(1).value;
            if (!genomeRaw || typeof genomeRaw !== 'string') return;
            if (!GENOME_ALIAS[genomeRaw.toLowerCase().trim()]) return;
            const chromRaw = row.getCell(3).value;
            const startRaw = row.getCell(4).value;
            const endRaw   = row.getCell(5).value;
            const refRaw   = row.getCell(6).value;
            const altRaw   = row.getCell(7).value;
            const posStart = parseInt(startRaw, 10);
            const posEnd   = parseInt(endRaw, 10) || posStart;
            if (!chromRaw || isNaN(posStart)) return;
            dataRows.push({
                rowNum,
                genomeStr: genomeRaw.trim(),
                chromStr:  String(chromRaw).trim(),
                posStart,
                posEnd,
                ref: refRaw ? String(refRaw).trim() : null,
                alt: altRaw ? String(altRaw).trim() : null,
            });
        });

        if (dataRows.length === 0) {
            throw new Error('Aucune ligne reconnue. Vérifiez que la colonne A contient un génome connu (ex : CanFam3.1).');
        }

        const errors = [];
        const t0 = Date.now();

        for (let i = 0; i < dataRows.length; i++) {
            const { rowNum, genomeStr, chromStr, posStart, posEnd, ref, alt } = dataRows[i];

            // Estimation temps restant
            const elapsed = (Date.now() - t0) / 1000;
            const avgSec  = i > 0 ? elapsed / i : 0;
            const remain  = i > 0 ? Math.round(avgSec * (dataRows.length - i)) : '…';
            msg.textContent = `Traitement ${i + 1} / ${dataRows.length} (ligne ${rowNum})${i > 0 ? ` — ~${remain}s restantes` : ''}`;
            bar.style.width = `${Math.round(((i) / dataRows.length) * 100)}%`;

            try {
                const { tokens, sequence } = await batchComputeMasked(genomeStr, chromStr, posStart, posEnd);
                const maskedStr = tokens.join('');
                // Remplacer le label [?] par [REF/ALT] issu des colonnes F/G si disponible
                const label = (ref && alt) ? `[${ref}/${alt}]` : null;
                const finalStr = label
                    ? maskedStr.replace(/\[[^\]]+\]/, label)
                    : maskedStr;
                const row = ws.getRow(rowNum);
                row.getCell(11).value = buildBatchRichText(finalStr);
                // Colonne N (14) : % GC de la séquence brute
                const gcCount = (sequence.match(/[GC]/g) || []).length;
                row.getCell(14).value = parseFloat(((gcCount / sequence.length) * 100).toFixed(1));
            } catch (e) {
                const marker = ws.getRow(rowNum).getCell(2).value || `ligne ${rowNum}`;
                errors.push({ rowNum, marker, error: e.message });
            }

            await new Promise(r => setTimeout(r, 300));
        }

        _batchOutputBuffer = await workbook.xlsx.writeBuffer();

        bar.style.width = '100%';
        msg.textContent = `Terminé — ${dataRows.length - errors.length} / ${dataRows.length} lignes traitées avec succès.`;
        showEl('batchDownload');

        if (errors.length > 0) {
            const html = `<p style="font-size:.83rem;font-weight:600;color:#b91c1c;margin-bottom:.5rem;">${errors.length} erreur(s) :</p>` +
                errors.map(e =>
                    `<div class="batch-error-item"><strong>${e.marker}</strong> (ligne ${e.rowNum}) — ${e.error}</div>`
                ).join('');
            document.getElementById('batchErrors').innerHTML = html;
            showEl('batchErrors');
        }

    } catch (e) {
        bar.style.width = '100%';
        bar.style.background = 'var(--mut-color)';
        msg.textContent = '⚠ ' + e.message;
    } finally {
        document.getElementById('batchRunBtn').disabled = false;
    }
}

function downloadBatch() {
    if (!_batchOutputBuffer) return;
    const blob = new Blob([_batchOutputBuffer], {
        type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'VariantScope_batch_output.xlsx';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 5000);
}

// Construit un objet richText ExcelJS : texte noir + [REF/VAR] en rouge
function buildBatchRichText(fullStr) {
    const match = fullStr.match(/\[[^\]]+\]/);
    if (!match) return { richText: [{ text: fullStr, font: { color: { argb: 'FF000000' } } }] };
    const before  = fullStr.slice(0, match.index);
    const bracket = match[0];
    const after   = fullStr.slice(match.index + bracket.length);
    const richText = [];
    if (before)  richText.push({ text: before,  font: { color: { argb: 'FF000000' } } });
    richText.push(              { text: bracket, font: { color: { argb: 'FFFF0000' } } });
    if (after)   richText.push({ text: after,   font: { color: { argb: 'FF000000' } } });
    return { richText };
}

// Point d'entrée batch : résout le génome et délègue au pipeline approprié
// Retourne { tokens, sequence }
async function batchComputeMasked(genomeStr, chromStr, posStart, posEnd) {
    const key     = genomeStr.toLowerCase().trim();
    const mapping = GENOME_ALIAS[key];
    if (!mapping) throw new Error(`Génome non reconnu : "${genomeStr}"`);

    const species   = SPECIES_CONFIG[mapping.species];
    const genomeCfg = species.genomes.find(g => g.id === mapping.genomeId);
    if (!genomeCfg) throw new Error(`Configuration introuvable pour ${mapping.genomeId}`);

    if (genomeCfg.backend === 'ucsc') {
        return await batchUCSCMasked(species, genomeCfg, chromStr, posStart, posEnd, 200);
    } else {
        return await batchEnsemblMasked(species, genomeCfg, chromStr, posStart, posEnd, 200);
    }
}

// Pipeline UCSC sans DOM — même logique que runUCSCPipeline, retourne tokens[]
async function batchUCSCMasked(species, genomeCfg, chrom, posStart, posEnd, windowSize) {
    // Normaliser le préfixe chr (CHR34 → chr34, 34 → chr34)
    chrom = 'chr' + chrom.replace(/^chr/i, '');

    let finalDb       = genomeCfg.ucscDb;
    let finalChrom    = chrom;
    let finalPosStart = posStart;
    let finalPosEnd   = posEnd;

    if (genomeCfg.ucscTarget && genomeCfg.ucscDb !== genomeCfg.ucscTarget) {
        const lifted  = await ucscLiftOver(genomeCfg.ucscDb, genomeCfg.ucscTarget, chrom, posStart);
        finalDb       = genomeCfg.ucscTarget;
        finalChrom    = lifted.chrom;
        finalPosStart = lifted.position;
        finalPosEnd   = lifted.position + (posEnd - posStart);
    }

    const baseChrom = resolveChromName(finalDb, finalChrom);
    const ucscChrom = /^NC_/.test(finalChrom) ? 'chr' + baseChrom : finalChrom;

    const midPos      = Math.round((finalPosStart + finalPosEnd) / 2);
    const winStart0   = Math.max(0, midPos - 1 - windowSize);
    const winEnd0     = midPos + windowSize;
    const winStart1   = winStart0 + 1;
    const mutIdxStart = finalPosStart - 1 - winStart0;
    const mutIdxEnd   = finalPosEnd   - 1 - winStart0;

    const sequence = await ucscSequence(finalDb, ucscChrom, winStart0, winEnd0);

    let variants = [];
    const GCF_ENSEMBL_MAP = { 'GCF_014441545.1': 'ROS_Cfam_1.0' };
    if (finalDb in GCF_ENSEMBL_MAP) {
        const ensChrom = baseChrom.replace(/^chr/i, '');
        variants = await ensemblVariants(species.ensemblSpecies, ensChrom, winStart1, winEnd0);
    } else if (genomeCfg.variantTrack) {
        variants = await ucscVariants(finalDb, genomeCfg.variantTrack, ucscChrom, winStart0, winEnd0);
    }

    const { tokens } = buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
    return { tokens, sequence };
}

// Pipeline Ensembl sans DOM — même logique que runEnsemblPipeline, retourne tokens[]
async function batchEnsemblMasked(species, genomeCfg, chrom, posStart, posEnd, windowSize) {
    let chromEns    = chrom.replace(/^chr/i, '');
    let finalChrom  = chromEns;
    let finalPosStart = posStart;
    let finalPosEnd   = posEnd;

    if (genomeCfg.ensemblAsm !== species.targetEnsemblAsm) {
        const lifted = await ensemblLiftOver(
            species.ensemblSpecies, genomeCfg.ensemblAsm,
            species.targetEnsemblAsm, chromEns, posStart
        );
        finalChrom    = lifted.chrom;
        finalPosStart = lifted.position;
        finalPosEnd   = lifted.position + (posEnd - posStart);
    }

    const midPos      = Math.round((finalPosStart + finalPosEnd) / 2);
    const winStart1   = Math.max(1, midPos - windowSize);
    const winEnd1     = midPos + windowSize;
    const mutIdxStart = finalPosStart - winStart1;
    const mutIdxEnd   = finalPosEnd   - winStart1;

    const sequence = await ensemblSequence(species.ensemblSpecies, finalChrom, winStart1, winEnd1);
    const variants = await ensemblVariants(species.ensemblSpecies, finalChrom, winStart1, winEnd1);

    const { tokens } = buildMaskedSeq(sequence, variants, winStart1, mutIdxStart, mutIdxEnd);
    return { tokens, sequence };
}

// =============================================================
//  Démarrage
// =============================================================
initForm();
