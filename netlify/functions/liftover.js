const https = require('https');

// Fait un GET HTTPS en suivant les redirections HTTPS uniquement
function httpsGet(url, maxRedirects = 5) {
    return new Promise((resolve, reject) => {
        https.get(url, (res) => {
            const { statusCode, headers } = res;

            if ([301, 302, 303, 307, 308].includes(statusCode)) {
                const location = headers.location || '';
                if (!location.startsWith('https://')) {
                    // Redirection vers HTTP ou page d'aide → base invalide
                    return resolve({ statusCode: 404, body: JSON.stringify({
                        error: `Assemblage non disponible dans l'API UCSC (redirigé vers ${location})`
                    })});
                }
                if (maxRedirects === 0) {
                    return resolve({ statusCode: 500, body: JSON.stringify({ error: 'Trop de redirections' }) });
                }
                return httpsGet(location, maxRedirects - 1).then(resolve).catch(reject);
            }

            let data = '';
            res.on('data', chunk => data += chunk);
            res.on('end', () => resolve({ statusCode, body: data }));
        }).on('error', (e) => reject(e));
    });
}

exports.handler = async (event) => {
    try {
        const params = new URLSearchParams(event.queryStringParameters).toString();
        const url = `https://api.genome.ucsc.edu/liftover?${params}`;
        const result = await httpsGet(url);

        return {
            statusCode: result.statusCode,
            headers: {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
            },
            body: result.body
        };
    } catch (e) {
        return {
            statusCode: 500,
            headers: { 'Content-Type': 'application/json', 'Access-Control-Allow-Origin': '*' },
            body: JSON.stringify({ error: e.message })
        };
    }
};
