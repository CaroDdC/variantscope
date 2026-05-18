const https = require('https');

module.exports = async function(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Content-Type', 'application/json');

    https.get('https://api.genome.ucsc.edu/list/ucscGenomes', (r) => {
        let data = '';
        r.on('data', chunk => data += chunk);
        r.on('end', () => res.status(200).send(data));
    }).on('error', (e) => res.status(500).json({ error: e.message }));
};
