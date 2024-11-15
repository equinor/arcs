import { useIsAuthenticated } from "@azure/msal-react";
import { useMsal } from "@azure/msal-react";
import { loginRequest } from "../api/AuthConfig";
import { Button } from "@equinor/eds-core-react";
import { IPublicClientApplication } from "@azure/msal-browser";

async function handleLogin(instance: IPublicClientApplication) {
    await instance.handleRedirectPromise();
    const accounts = instance.getAllAccounts();
    if (accounts.length === 0) {
        await instance.loginRedirect(loginRequest).catch((e: any) => {
            console.error(e);
        });
    }
}

export const SignInButton = () => {
    const { instance } = useMsal();

    return (
        <Button href="" variant="contained" onClick={() => handleLogin(instance)}>
            Sign in using Redirect
        </Button>
    );
};

export const SignInPage = () => {
    const isAuthenticated = useIsAuthenticated();

    return (
        <>
            <div className="sign-in-button">{isAuthenticated ? <span>Signed In</span> : <SignInButton />}</div>
        </>
    );
};