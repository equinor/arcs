import { Configuration } from "@azure/msal-browser";
import { IPublicClientApplication } from "@azure/msal-browser";
import config from "../configuration";

export const msalConfig: Configuration = {
    auth: {
        clientId: config.AD_CLIENT_ID,
        authority: "https://login.microsoftonline.com/" + config.AD_TENANT_ID,
        redirectUri: window.location.origin,
    },
    cache: {
        cacheLocation: "sessionStorage",
    },
};

export const loginRequest = {
    scopes: [config.BACKEND_API_SCOPE],
};

export async function fetchAccessToken(
    instance: IPublicClientApplication,
    environment: string,
    homeAccountId: string,
    username: string,
    tenantId: string,
    localAccountId: string
): Promise<string> {
    // Silently acquires an access token which is then attached to a request for Microsoft Graph data
    return instance
        .acquireTokenSilent({
            ...loginRequest,
            account: {
                environment: environment,
                homeAccountId: homeAccountId,
                username: username,
                tenantId: tenantId,
                localAccountId: localAccountId,
            },
        })
        .then((response: any) => {
            const accessToken: string = response.accessToken ?? "";
            return accessToken;
        })
        .catch((e: any) => {
            console.log(e);
            return instance.acquireTokenRedirect(loginRequest).then(() => {
                return "The page should be refreshed automatically.";
            });
        });
}
