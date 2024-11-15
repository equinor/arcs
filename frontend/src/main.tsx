/* eslint-disable @typescript-eslint/no-explicit-any */
import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App";
import { PublicClientApplication } from "@azure/msal-browser";
import { MsalProvider } from "@azure/msal-react";
import { msalConfig } from "./api/AuthConfig";
import store from "./store/Store";
import { Provider } from "react-redux";
// ========================================

export const msalInstance = new PublicClientApplication(msalConfig);

// Error handling for redirect login
msalInstance.initialize().then(() => {
    msalInstance.handleRedirectPromise().catch((error: any) => {
        console.error(error);
    });
}); 

ReactDOM.createRoot(document.getElementById("root") as HTMLElement).render(
    <React.StrictMode>
        <MsalProvider instance={msalInstance}>
            <Provider store={store}>
                <App />
            </Provider>
        </MsalProvider>
    </React.StrictMode>
);
