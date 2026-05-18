// Liftover via fichiers chain UCSC (streamés et parsés à la volée)
// Évite les limites de l'API REST UCSC qui ne connaît pas tous les assemblages

const https = require('https');
const zlib  = require('zlib');

// ----------- Téléchargement streamé du fichier chain --------

function liftoverWithChainFile(fromDb, toDb, chrom, pos0) {
    // Convention UCSC pour les noms de fichiers chain
    const toDbCap = toDb.charAt(0).toUpperCase() + toDb.slice(1);
    const url = `https://hgdownload.soe.ucsc.edu/goldenPath/${fromDb}/liftOver/` +
                `${fromDb}To${toDbCap}.over.chain.gz`;

    return new Promise((resolve, reject) => {
        https.get(url, (res) => {
            if (res.statusCode === 404) {
                return reject(new Error(
                    `Fichier chain introuvable : ${fromDb}→${toDb}. ` +
                    `URL essayée : ${url}`
                ));
            }
            if (res.statusCode !== 200) {
                return reject(new Error(`Téléchargement chain : HTTP ${res.statusCode}`));
            }

            const gunzip = zlib.createGunzip();
            res.pipe(gunzip);

            let buffer       = '';
            let chain        = null;  // chaîne courante
            let tPos, qPos;
            let result       = null;
            let done         = false;

            function processLine(line) {
                if (done) return;
                const t = line.trim();
                if (!t) return;

                if (t.startsWith('chain')) {
                    // En-tête d'une nouvelle chaîne
                    const p = t.split(/\s+/);
                    chain = {
                        tName:   p[2],  // chrom source
                        tStrand: p[4],
                        tStart:  parseInt(p[5]),
                        tEnd:    parseInt(p[6]),
                        qName:   p[7],  // chrom cible
                        qSize:   parseInt(p[8]),
                        qStrand: p[9],
                        qStart:  parseInt(p[10])
                    };
                    tPos = chain.tStart;
                    qPos = chain.qStart;
                } else if (chain) {
                    // Bloc d'alignement
                    // On ne traite que si le chrom source correspond et couvre notre position
                    if (chain.tName !== chrom || pos0 < chain.tStart || pos0 >= chain.tEnd) return;

                    const p    = t.split(/\s+/);
                    const size = parseInt(p[0]);
                    const dt   = p[1] ? parseInt(p[1]) : 0;
                    const dq   = p[2] ? parseInt(p[2]) : 0;

                    if (pos0 >= tPos && pos0 < tPos + size) {
                        const offset = pos0 - tPos;
                        let mappedPos;
                        if (chain.qStrand === '+') {
                            mappedPos = qPos + offset;
                        } else {
                            // Brin complémentaire : coordonnée depuis la fin
                            mappedPos = chain.qSize - (qPos + size - offset) - 1;
                        }
                        result = { chrom: chain.qName, pos0: mappedPos, strand: chain.qStrand };
                        done = true;
                        try { res.destroy(); } catch {}  // Arrêter le téléchargement
                    } else {
                        tPos += size + dt;
                        qPos += size + dq;
                    }
                }
            }

            gunzip.on('data', (chunk) => {
                if (done) return;
                buffer += chunk.toString('ascii');
                const lines = buffer.split('\n');
                buffer = lines.pop();  // Garder la ligne incomplète
                lines.forEach(processLine);
            });

            const finish = () => resolve(result);
            gunzip.on('end',   finish);
            gunzip.on('error', (e) => done ? resolve(result) : reject(e));
            res.on('error',    (e) => done ? resolve(result) : reject(e));
            res.on('close',    finish);
        }).on('error', reject);
    });
}

// ----------- Handler Netlify --------------------------------

exports.handler = async (event) => {
    const headers = {
        'Content-Type': 'application/json',
        'Access-Control-Allow-Origin': '*'
    };

    const { fromDb, toDb, chrom, start, end } = event.queryStringParameters || {};
    if (!fromDb || !toDb || !chrom || start === undefined) {
        return { statusCode: 400, headers,
                 body: JSON.stringify({ error: 'Paramètres manquants (fromDb, toDb, chrom, start, end)' }) };
    }

    const pos0 = parseInt(start);

    try {
        const result = await liftoverWithChainFile(fromDb, toDb, chrom, pos0);

        if (!result) {
            return { statusCode: 200, headers, body: JSON.stringify({
                mappedCoordinates:   [],
                unmappedCoordinates: [{ reason: 'Position non mappable dans cet assemblage cible' }]
            })};
        }

        // Format identique à l'API REST UCSC
        return { statusCode: 200, headers, body: JSON.stringify({
            mappedCoordinates: [{
                chrom:  result.chrom,
                start:  result.pos0,
                end:    result.pos0 + 1,
                strand: result.strand
            }],
            unmappedCoordinates: []
        })};

    } catch (e) {
        return { statusCode: 500, headers, body: JSON.stringify({ error: e.message }) };
    }
};
