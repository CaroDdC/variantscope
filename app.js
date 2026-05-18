// =============================================================
//  VariantScope — app.js
//  Dépendances : aucune (vanilla JS)
//  APIs utilisées :
//    - UCSC REST API  : séquence + LiftOver
//    - Ensembl REST API : polymorphismes
// =============================================================

// ----------- Configuration des espèces / génomes -------------

const SPECIES_CONFIG = {
    human: {
        label:          'Humain (Homo sapiens)',
        ensemblSpecies: 'homo_sapiens',
        genomes: [
            { id: 'hg38',    label: 'hg38 / GRCh38 (récent)' },
            { id: 'hg19',    label: 'hg19 / GRCh37' }
        ],
        // Génome cible pour Ensembl (doit correspondre à l'assembly Ensembl par défaut)
        targetGenome: 'hg38'
    },
    mouse: {
        label:          'Souris (Mus musculus)',
        ensemblSpecies: 'mus_musculus',
        genomes: [
            { id: 'mm39',    label: 'mm39 / GRCm39 (récent)' },
            { id: 'mm10',    label: 'mm10 / GRCm38' }
        ],
        targetGenome: 'mm39'
    },
    dog: {
        label:          'Chien (Canis lupus familiaris)',
        ensemblSpecies: 'canis_lupus_familiaris',
        genomes: [
            { id: 'canFam5', label: 'canFam5 / ROS_Cfam_1.0 (récent)' },
            { id: 'canFam3', label: 'canFam3 / CanFam3.1' }
        ],
        // canFam5 = ROS_Cfam_1.0 : assembly utilisé par Ensembl
        targetGenome: 'canFam5'
    },
    cat: {
        label:          'Chat (Felis catus)',
        ensemblSpecies: 'felis_catus',
        genomes: [
            { id: 'felCat9', label: 'felCat9 / Felis_catus_9.0 (récent)' },
            { id: 'felCat8', label: 'felCat8 / Felis_catus_8.0' }
        ],
        // felCat9 = Felis_catus_9.0 : assembly utilisé par Ensembl
        targetGenome: 'felCat9'
    },
    horse: {
        label:          'Cheval (Equus caballus)',
        ensemblSpecies: 'equus_caballus',
        genomes: [
            { id: 'equCab3', label: 'equCab3 / EquCab3.0 (récent)' },
            { id: 'equCab2', label: 'equCab2 / EquCab2.0' }
        ],
        // equCab3 = EquCab3.0 : assembly utilisé par Ensembl
        targetGenome: 'equCab3'
    }
};

// Séquences brutes stockées pour le bouton "Copier"
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
    const genome     = document.getElementById('genome').value;
    let   chrom      = document.getElementById('chromosome').value.trim();
    const position   = parseInt(document.getElementById('position').value, 10);
    const windowSize = Math.max(50, parseInt(document.getElementById('windowSize').value, 10) || 500);

    // Normaliser le préfixe "chr"
    if (!/^chr/i.test(chrom)) chrom = 'chr' + chrom;

    const btn = document.getElementById('submitBtn');
    btn.disabled = true;
    resetUI();
    showEl('results');
    showEl('loading');
    hideEl('errorBox');

    try {
        const species = SPECIES_CONFIG[speciesKey];
        let finalGenome = genome;
        let finalChrom  = chrom;
        let finalPos    = position;

        // --- Étape 1 : LiftOver (si le génome d'entrée ≠ génome cible Ensembl) ---
        if (genome !== species.targetGenome) {
            setLoadingMsg('LiftOver en cours…');
            const lifted = await doLiftOver(genome, species.targetGenome, chrom, position);
            finalGenome = species.targetGenome;
            finalChrom  = lifted.chrom;
            finalPos    = lifted.position;
            showLiftoverBanner(genome, species.targetGenome, chrom, position, finalChrom, finalPos);
        }

        // --- Étape 2 : Coordonnées de la fenêtre (UCSC = 0-based half-open) ---
        const winStart0 = Math.max(0, finalPos - 1 - windowSize); // inclusif, 0-based
        const winEnd0   = finalPos + windowSize;                   // exclusif, 0-based

        // --- Étape 3 : Séquence ---
        setLoadingMsg('Récupération de la séquence UCSC…');
        const sequence = await fetchSequence(finalGenome, finalChrom, winStart0, winEnd0);

        // --- Étape 4 : Variants Ensembl ---
        setLoadingMsg('Récupération des variants Ensembl…');
        const chromEns  = finalChrom.replace(/^chr/i, '');  // Ensembl n'utilise pas "chr"
        const regionStart1 = winStart0 + 1;                 // passage en 1-based
        const regionEnd1   = winEnd0;                       // winEnd0 est exclusif → correspond à la dernière base en 1-based
        const variants = await fetchVariants(species.ensemblSpecies, chromEns, regionStart1, regionEnd1);

        // --- Étape 5 : Rendu ---
        // Index 0-based de la position de mutation dans la chaîne `sequence`
        const mutIdx = finalPos - 1 - winStart0;

        _rawSeq = sequence;
        displaySequence(sequence, mutIdx, finalChrom, regionStart1, regionEnd1);
        displayVariants(variants);

        const masked = buildMaskedSeq(sequence, variants, regionStart1);
        _maskedSeq = masked;
        displayMaskedSequence(masked, variants, mutIdx, finalChrom, regionStart1, regionEnd1);

    } catch (err) {
        showError(err.message || 'Erreur inattendue. Vérifiez vos paramètres.');
    } finally {
        hideEl('loading');
        btn.disabled = false;
    }
});

// ----------- Appels API --------------------------------------

async function doLiftOver(fromDb, toDb, chrom, position) {
    // UCSC LiftOver API — coordonnées 0-based
    const start = position - 1;
    const end   = position;
    const url = `https://api.genome.ucsc.edu/liftover` +
                `?fromDb=${encodeURIComponent(fromDb)}` +
                `&toDb=${encodeURIComponent(toDb)}` +
                `&chrom=${encodeURIComponent(chrom)}` +
                `&start=${start}&end=${end}`;

    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`LiftOver : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    if (data.error) throw new Error('LiftOver : ' + data.error);

    if (!data.mappedCoordinates || data.mappedCoordinates.length === 0) {
        const reason = data.unmappedCoordinates?.[0]?.reason || 'position non mappable dans ce génome cible';
        throw new Error(`LiftOver : aucune correspondance trouvée (${reason})`);
    }

    const m = data.mappedCoordinates[0];
    return {
        chrom:    m.chrom,
        position: m.start + 1  // retour en 1-based
    };
}

async function fetchSequence(genome, chrom, start0, end0) {
    const url = `https://api.genome.ucsc.edu/getData/sequence` +
                `?genome=${encodeURIComponent(genome)}` +
                `&chrom=${encodeURIComponent(chrom)}` +
                `&start=${start0}&end=${end0}`;

    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`Séquence : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    if (data.error) throw new Error('Séquence UCSC : ' + data.error);

    return data.dna.toUpperCase();
}

async function fetchVariants(ensemblSpecies, chromEns, start1, end1) {
    const url = `https://rest.ensembl.org/overlap/region/${ensemblSpecies}` +
                `/${chromEns}:${start1}-${end1}` +
                `?feature=variation;content-type=application/json`;

    const resp = await fetch(url);
    if (resp.status === 404) return [];  // région sans variants connus
    if (!resp.ok) throw new Error(`Variants Ensembl : erreur HTTP ${resp.status}`);

    const data = await resp.json();
    return Array.isArray(data) ? data : [];
}

// ----------- Construction de la séquence masquée -------------

function buildMaskedSeq(seq, variants, regionStart1) {
    // regionStart1 : position 1-based de la première base de `seq`
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
        const shown = variants.slice(0, 500);
        const rows  = shown.map(v => {
            const name   = v.id || v.variation_name || '–';
            const link   = /^rs\d+/.test(name)
                           ? `<a href="https://www.ensembl.org/id/${name}" target="_blank" rel="noopener">${name}</a>`
                           : name;
            const alleles = Array.isArray(v.alleles) ? v.alleles.join('/') :
                            (typeof v.alleles === 'string' ? v.alleles : '–');
            const conseq  = Array.isArray(v.consequence_type)
                            ? v.consequence_type[0]
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

    // Construire la map des highlights
    const hMap = new Map();
    for (let i = 0; i < masked.length; i++) {
        if (masked[i] === 'N') hMap.set(i, 'n');
    }
    hMap.set(mutIdx, 'pos'); // la position de mutation prend priorité

    document.getElementById('maskedDisplay').innerHTML = renderSeq(masked, hMap);
    showEl('maskedCard');
}

// ----------- Rendu HTML d'une séquence -----------------------

function renderSeq(seq, highlightMap) {
    const LINE  = 60;  // bases par ligne
    const BLOCK = 10;  // espace tous les N bases
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
        // Feedback visuel temporaire
        document.querySelectorAll('.copy-btn').forEach(btn => {
            if (btn.getAttribute('onclick')?.includes(which)) {
                btn.textContent = 'Copié !';
                setTimeout(() => { btn.textContent = 'Copier'; }, 1500);
            }
        });
    }).catch(() => {
        alert('La copie automatique est bloquée par votre navigateur.\nSélectionnez la séquence manuellement.');
    });
}

// ----------- Helpers d'interface ----------------------------

function showEl(id) { document.getElementById(id)?.classList.remove('hidden'); }
function hideEl(id) { document.getElementById(id)?.classList.add('hidden'); }

function resetUI() {
    ['sequenceCard', 'variantsCard', 'maskedCard', 'liftoverInfo', 'loading']
        .forEach(hideEl);
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

function showLiftoverBanner(fromDb, toDb, fromChrom, fromPos, toChrom, toPos) {
    const el = document.getElementById('liftoverInfo');
    el.innerHTML =
        `<strong>LiftOver effectué :</strong> ` +
        `${fromDb} ${fromChrom}:${fmt(fromPos)} ` +
        `→ ${toDb} ${toChrom}:${fmt(toPos)}`;
    showEl('liftoverInfo');
}

function fmt(n) {
    return typeof n === 'number' ? n.toLocaleString('fr-FR') : n;
}

// ----------- Démarrage --------------------------------------

initForm();
