const https = require('https');

exports.handler = async (event) => {
    const params = new URLSearchParams(event.queryStringParameters).toString();
    const url = `https://api.genome.ucsc.edu/liftover?${params}`;

    return new Promise((resolve) => {
        https.get(url, (res) => {
            let data = '';
            res.on('data', chunk => data += chunk);
            res.on('end', () => resolve({
                statusCode: res.statusCode,
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
