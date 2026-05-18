const https = require('https');
const zlib  = require('zlib');

function liftoverWithChainFile(fromDb, toDb, chrom, pos0) {
    const toDbCap = toDb.charAt(0).toUpperCase() + toDb.slice(1);
    const url = `https://hgdownload.soe.ucsc.edu/goldenPath/${fromDb}/liftOver/` +
                `${fromDb}To${toDbCap}.over.chain.gz`;

    return new Promise((resolve, reject) => {
        https.get(url, (res) => {
            if (res.statusCode === 404) {
                return reject(new Error(`Fichier chain introuvable : ${fromDb}→${toDb}. URL essayée : ${url}`));
            }
            if (res.statusCode !== 200) {
                return reject(new Error(`Téléchargement chain : HTTP ${res.statusCode}`));
            }

            const gunzip = zlib.createGunzip();
            res.pipe(gunzip);

            let buffer = '';
            let chain  = null;
            let tPos, qPos;
            let result = null;
            let done   = false;

            function processLine(line) {
                if (done) return;
                const t = line.trim();
                if (!t) return;

                if (t.startsWith('chain')) {
                    const p = t.split(/\s+/);
                    chain = {
                        tName: p[2], tStrand: p[4],
                        tStart: parseInt(p[5]), tEnd: parseInt(p[6]),
                        qName: p[7], qSize: parseInt(p[8]),
                        qStrand: p[9], qStart: parseInt(p[10])
                    };
                    tPos = chain.tStart;
                    qPos = chain.qStart;
                } else if (chain) {
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
                            mappedPos = chain.qSize - (qPos + size - offset) - 1;
                        }
                        result = { chrom: chain.qName, pos0: mappedPos, strand: chain.qStrand };
                        done = true;
                        try { res.destroy(); } catch {}
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
                buffer = lines.pop();
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

module.exports = async function(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Content-Type', 'application/json');

    const { fromDb, toDb, chrom, start } = req.query;
    if (!fromDb || !toDb || !chrom || start === undefined) {
        return res.status(400).json({ error: 'Paramètres manquants (fromDb, toDb, chrom, start)' });
    }

    const pos0 = parseInt(start);

    try {
        const result = await liftoverWithChainFile(fromDb, toDb, chrom, pos0);

        if (!result) {
            return res.status(200).json({
                mappedCoordinates:   [],
                unmappedCoordinates: [{ reason: 'Position non mappable dans cet assemblage cible' }]
            });
        }

        return res.status(200).json({
            mappedCoordinates: [{
                chrom:  result.chrom,
                start:  result.pos0,
                end:    result.pos0 + 1,
                strand: result.strand
            }],
            unmappedCoordinates: []
        });
    } catch (e) {
        return res.status(500).json({ error: e.message });
    }
};
