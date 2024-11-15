interface EsiTronConfiguration {
    BACKEND_URL: string;
    BACKEND_API_SCOPE: string;
    AD_CLIENT_ID: string;
    AD_TENANT_ID: string;
    FEEDBACK_FORM_URL: string;
    APPLICATION_INSIGHTS_CONNECTION_STRING: string;
}

declare global {
    interface Window {
        injectEnv: EsiTronConfiguration;
    }

    interface ImportMetaEnv {
        DEV: boolean;
        VITE_BACKEND_URL: string;
        VITE_BACKEND_API_SCOPE: string;
        VITE_AD_CLIENT_ID: string;
        VITE_AD_TENANT_ID: string;
        VITE_FEEDBACK_FORM_URL: string;
        VITE_APPLICATION_INSIGHTS_CONNECTION_STRING: string;
    }

    interface ImportMeta {
        readonly env: ImportMetaEnv;
    }
}

function getEnvVars(): EsiTronConfiguration {
    if (import.meta.env.DEV) {
        // Used for local development only
        const config: EsiTronConfiguration = {
            BACKEND_URL: import.meta.env.VITE_BACKEND_URL,
            BACKEND_API_SCOPE: import.meta.env.VITE_BACKEND_API_SCOPE,
            AD_CLIENT_ID: import.meta.env.VITE_AD_CLIENT_ID,
            AD_TENANT_ID: import.meta.env.VITE_AD_TENANT_ID,
            FEEDBACK_FORM_URL: import.meta.env.VITE_FEEDBACK_FORM_URL,
            APPLICATION_INSIGHTS_CONNECTION_STRING: import.meta.env.VITE_APPLICATION_INSIGHTS_CONNECTION_STRING,
        };
        return config;
    }

    // injectEnv is injected by inject-env.js based on inject-env-template.js
    return window.injectEnv;
}

const config: EsiTronConfiguration = getEnvVars();

export default config;