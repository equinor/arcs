import { AuthenticatedTemplate, UnauthenticatedTemplate } from "@azure/msal-react";
import { SignInPage } from "./pages/SignInPage";
import { AppInsightsContext } from "@microsoft/applicationinsights-react-js";
import { MainPage } from "./pages/MainPage";
import "./App.css";
import { reactPlugin } from "./utils/appplicationInsights";

function App() {
    return (
        <>
            <UnauthenticatedTemplate>
                <div className="sign-in-page">
                    <SignInPage></SignInPage>
                </div>
            </UnauthenticatedTemplate>
            <AuthenticatedTemplate>
                {/* <AppInsightsContext.Provider value={reactPlugin}> */}
                    <MainPage />
                {/* </AppInsightsContext.Provider> */}
            </AuthenticatedTemplate>
        </>
    );
}

export default App;

