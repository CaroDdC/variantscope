const https = require('https');

exports.handler = async () => {
    return new Promise((resolve) => {
        https.get('https://api.genome.ucsc.edu/list/ucscGenomes', (res) => {
            let data = '';
            res.on('data', chunk => data += chunk);
            res.on('end', () => resolve({
                statusCode: 200,
                headers: {
                    'Content-Type': 'application/json',
                    'Access-Control-Allow-Origin': '*'
                },
                body: data
            }));
        }).on('error', (e) => resolve({
            statusCode: 500,
            body: JSON.stringify({ error: e.message })
        }));
    });
};
