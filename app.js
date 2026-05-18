// =============================================================
//  VariantScope — app.js
//  API utilisée : Ensembl REST API uniquement (supporte CORS)
//    - /map/{species}/{asm1}/{region}/{asm2} → LiftOver
//    - /sequence/region/{species}/{region}   → Séquence
//    - /overlap/region/{species}/{region}    → Variants
// =============================================================

// ----------- Configuration des espèces / génomes -------------
// ensemblAsm : nom de l'assemblage tel qu'Ensembl le connaît
// targetEnsemblAsm : assemblage cible (celui indexé par Ensembl pour les variants)

const SPECIES_CONFIG = {
    dog: {
        label:            'Chien (Canis lupus familiaris)',
        ensemblSpecies:   'canis_lupus_familiaris',
        genomes: [
            { id: 'ROS_Cfam_1.0', label: 'ROS_Cfam_1.0 / canFam5 (cible Ensembl)', ensemblAsm: 'ROS_Cfam_1.0' },
            { id: 'CanFam3.1',    label: 'CanFam3.1 / canFam3',                      ensemblAsm: 'CanFam3.1'    }
        ],
        targetEnsemblAsm: 'ROS_Cfam_1.0'
    },
    cat: {
        label:            'Chat (Felis catus)',
        ensemblSpecies:   'felis_catus',
        genomes: [
            { id: 'Felis_catus_9.0', label: 'Felis_catus_9.0 / felCat9 (récent)', ensemblAsm: 'Felis_catus_9.0' },
            { id: 'Felis_catus_8.0', label: 'Felis_catus_8.0 / felCat8',          ensemblAsm: 'Felis_catus_8.0' }
        ],
        targetEnsemblAsm: 'Felis_catus_9.0'
    },
    horse: {
        label:            'Cheval (Equus caballus)',
        ensemblSpecies:   'equus_caballus',
        genomes: [
            { id: 'EquCab3.0', label: 'EquCab3.0 / equCab3 (récent)', ensemblAsm: 'EquCab3.0' },
            { id: 'EquCab2.0', label: 'EquCab2.0 / equCab2',          ensemblAsm: 'EquCab2.0' }
        ],
        targetEnsemblAsm: 'EquCab3.0'
    }
};

// Séquences stockées pour le bouton "Copier"
let _rawSeq    = '';
let _maskedSeq = '';

// ----------- Initialisation du formulaire --------------------

function initForm() {
    const speciesSel = document.getElementById('species');
    Object.entries(SPECIES_CONFIG).forEach(([key, cfg]) => {
        const opt = document.createElement('option');
        opt.value = key;
        opt.textContent = cfg.label;
        speciesSel.appendChild(opt);
    });
    refreshGenomeOptions();
}

function refreshGenomeOptions() {
    const key = document.getElementById('species').value;
    const genomeSel = document.getElementById('genome');
    genomeSel.innerHTML = '';
    SPECIES_CONFIG[key].genomes.forEach(g => {
        const opt = document.createElement('option');
        opt.value = g.id;
        opt.textContent = g.label;
        genomeSel.appendChild(opt);
    });
}

document.getElementById('species').addEventListener('change', refreshGenomeOptions);

// ----------- Soumission du formulaire ------------------------

document.getElementById('variantForm').addEventListener('submit', async (e) => {
    e.preventDefault();

    const speciesKey = document.getElementById('species').value;
    const genomeId   = document.getElementById('genome').value;
    let   chrom      = document.getElementById('chromosome').value.trim();
    const position   = parseInt(document.getElementById('position').value, 10);
    const windowSize = Math.max(50, parseInt(document.getElementById('windowSize').value, 10) || 500);

    // Ensembl n'utilise pas le préfixe "chr"
    const chromEns = chrom.replace(/^chr/i, '');

    const btn = document.getElementById('submitBtn');
    btn.disabled = true;
    resetUI();
    showEl('results');
    showEl('loading');
    hideEl('errorBox');

    try {
        const species = SPECIES_CONFIG[speciesKey];

        let finalAsm   = genomeId;
        let finalChrom = chromEns;
        let finalPos   = position;

        // --- Étape 1 : LiftOver via Ensembl (si assemblage ≠ cible) ---
        if (genomeId !== species.targetEnsemblAsm) {
            setLoadingMsg('LiftOver en cours (Ensembl)…');
            const lifted = await doLiftOver(
                species.ensemblSpecies,
                genomeId,
                species.targetEnsemblAsm,
                chromEns,
                position
            );
            finalAsm   = species.targetEnsemblAsm;
            finalChrom = lifted.chrom;
            finalPos   = lifted.position;
            showLiftoverBanner(genomeId, species.targetEnsemblAsm, chromEns, position, finalChrom, finalPos);
        }

        // --- Étape 2 : Coordonnées de la fenêtre (1-based, Ensembl) ---
        const winStart1 = Math.max(1, finalPos - windowSize);
        const winEnd1   = finalPos + windowSize;

        // --- Étape 3 : Séquence ---
        setLoadingMsg('Récupération de la séquence (Ensembl)…');
        const sequence = await fetchSequence(species.ensemblSpecies, finalChrom, winStart1, winEnd1);

        // --- Étape 4 : Variants ---
        setLoadingMsg('Récupération des variants (Ensembl)…');
        const variants = await fetchVariants(species.ensemblSpecies, finalChrom, winStart1, winEnd1);

        // Index 0-based de la position dans la séquence
        const mutIdx = finalPos - winStart1;

        // --- Étape 5 : Affichage ---
        _rawSeq = sequence;
        displaySequence(sequence, mutIdx, finalChrom, winStart1, winEnd1);
        displayVariants(variants);

        const masked = buildMaskedSeq(sequence, variants, winStart1);
        _maskedSeq = masked;
        displayMaskedSequence(masked, variants, mutIdx, finalChrom, winStart1, winEnd1);

    } catch (err) {
        showError(err.message || 'Erreur inattendue. Vérifiez vos paramètres.');
    } finally {
        hideEl('loading');
        btn.disabled = false;
    }
});

// ----------- Appels API Ensembl (CORS OK) --------------------

async function doLiftOver(ensemblSpecies, fromAsm, toAsm, chrom, position) {
    // Ensembl /map : coordonnées 1-based inclusives
    const url = `https://rest.ensembl.org/map/${ensemblSpecies}` +
                `/${encodeURIComponent(fromAsm)}` +
                `/${encodeURIComponent(chrom)}:${position}..${position}` +
                `/${encodeURIComponent(toAsm)}` +
                `?content-type=application/json`;

    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`LiftOver : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    if (data.error) throw new Error('LiftOver : ' + data.error);
    if (!data.mappings || data.mappings.length === 0) {
        throw new Error('LiftOver : aucune correspondance trouvée pour cette position dans le génome cible.');
    }

    const m = data.mappings[0].mapped;
    return {
        chrom:    m.seq_region_name,   // sans préfixe "chr"
        position: m.start              // 1-based
    };
}

async function fetchSequence(ensemblSpecies, chrom, start1, end1) {
    // Ensembl /sequence/region : coordonnées 1-based inclusives
    const url = `https://rest.ensembl.org/sequence/region/${ensemblSpecies}` +
                `/${encodeURIComponent(chrom)}:${start1}..${end1}` +
                `?content-type=application/json`;

    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Séquence : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    if (data.error) throw new Error('Séquence : ' + data.error);

    return data.seq.toUpperCase();
}

async function fetchVariants(ensemblSpecies, chrom, start1, end1) {
    const url = `https://rest.ensembl.org/overlap/region/${ensemblSpecies}` +
                `/${encodeURIComponent(chrom)}:${start1}-${end1}` +
                `?feature=variation;content-type=application/json`;

    const resp = await fetch(url);
    if (resp.status === 404) return [];
    if (!resp.ok) throw new Error(`Variants : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    return Array.isArray(data) ? data : [];
}

// ----------- Construction de la séquence masquée -------------

function buildMaskedSeq(seq, variants, regionStart1) {
    const bases = seq.split('');
    variants.forEach(v => {
        const idxStart = v.start - regionStart1;
        const idxEnd   = (v.end !== undefined ? v.end : v.start) - regionStart1;
        for (let i = idxStart; i <= idxEnd; i++) {
            if (i >= 0 && i < bases.length) bases[i] = 'N';
        }
    });
    return bases.join('');
}

// ----------- Affichage des résultats -------------------------

function displaySequence(seq, mutIdx, chrom, start1, end1) {
    document.getElementById('coordsDisplay').textContent =
        `${chrom}:${fmt(start1)}–${fmt(end1)}  ·  ${seq.length} pb`;
    document.getElementById('sequenceDisplay').innerHTML =
        renderSeq(seq, new Map([[mutIdx, 'pos']]));
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
            const name   = v.id || v.variation_name || '–';
            const link   = /^rs\d+/.test(name)
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

function displayMaskedSequence(masked, variants, mutIdx, chrom, start1, end1) {
    document.getElementById('maskedCoordsDisplay').textContent =
        `${chrom}:${fmt(start1)}–${fmt(end1)}  ·  ${variants.length} position(s) masquée(s)`;

    const hMap = new Map();
    for (let i = 0; i < masked.length; i++) {
        if (masked[i] === 'N') hMap.set(i, 'n');
    }
    hMap.set(mutIdx, 'pos');

    document.getElementById('maskedDisplay').innerHTML = renderSeq(masked, hMap);
    showEl('maskedCard');
}

// ----------- Rendu HTML d'une séquence -----------------------

function renderSeq(seq, highlightMap) {
    const LINE  = 60;
    const BLOCK = 10;
    let html = '<div class="seq-display">';

    for (let ls = 0; ls < seq.length; ls += LINE) {
        const le = Math.min(ls + LINE, seq.length);
        html += `<div class="seq-line">`;
        html += `<span class="seq-coord">${ls + 1}</span>`;
        html += `<span class="seq-bases">`;

        for (let i = ls; i < le; i++) {
            if (i > ls && (i - ls) % BLOCK === 0) html += ' ';
            const base = seq[i];
            const type = highlightMap.get(i);
            if      (type === 'pos') html += `<span class="highlight-pos">${base}</span>`;
            else if (type === 'n')   html += `<span class="highlight-n">N</span>`;
            else                     html += base;
        }

        html += `</span></div>`;
    }
    return html + '</div>';
}

// ----------- Copie dans le presse-papiers -------------------

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
    }).catch(() => {
        alert('Copie automatique bloquée. Sélectionnez la séquence manuellement.');
    });
}

// ----------- Helpers interface -------------------------------

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

function showLiftoverBanner(fromAsm, toAsm, fromChrom, fromPos, toChrom, toPos) {
    const el = document.getElementById('liftoverInfo');
    el.innerHTML =
        `<strong>LiftOver effectué :</strong> ` +
        `${fromAsm} ${fromChrom}:${fmt(fromPos)} → ${toAsm} ${toChrom}:${fmt(toPos)}`;
    showEl('liftoverInfo');
}

function fmt(n) {
    return typeof n === 'number' ? n.toLocaleString('fr-FR') : n;
}

// ----------- Démarrage --------------------------------------
initForm();
