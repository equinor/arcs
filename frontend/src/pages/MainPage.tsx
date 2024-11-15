import { useMsal } from "@azure/msal-react";
import { fetchAccessToken } from "../api/AuthConfig";
import { Header } from "../components/header/Header";
import { createContext, useEffect } from "react";
import { BrowserRouter, Route, Routes } from "react-router-dom";
import styled from "styled-components";
import { Typography } from "@equinor/eds-core-react";
import ErrorDialog from "./ErrorDialog";
import { setCredentials } from "../store/AuthSlice";
import { useAppDispatch, useAppSelector } from "../store/Hooks";

const StyledSite = styled.div`
    display: flex;
    flex-direction: column;
    gap: 16px;
`;

export const AccessTokenContext = createContext(["", ""]);


export function MainPage() {
    const authContext = useMsal();
    const accessToken = useAppSelector((state) => state.auth.token);
    const user = useAppSelector((state) => state.auth.username);
    const dispatch = useAppDispatch();

    useEffect(() => {
        fetchAccessToken(
            authContext.instance,
            authContext.accounts[0].environment,
            authContext.accounts[0].homeAccountId,
            authContext.accounts[0].username,
            authContext.accounts[0].tenantId,
            authContext.accounts[0].localAccountId
        ).then((result) => {
            dispatch(
                setCredentials({
                    username: authContext.accounts[0].username,
                    token: result,
                    environment: authContext.accounts[0].environment,
                    homeAccountId: authContext.accounts[0].homeAccountId,
                    tenantId: authContext.accounts[0].tenantId,
                    localAccountId: authContext.accounts[0].localAccountId,
                })
            );
        });
    }, []);

    return (
        <>
            {accessToken === "" && <Typography>Loading...</Typography>}
            {accessToken !== "" && (
                <>
                    <AccessTokenContext.Provider value={[accessToken, user]}>
                        <BrowserRouter>
                            <StyledSite>
                                <Header />
                                <Routes>
                                    <Route path="/" element={<Typography>Welcome to the Co2Spec site</Typography>} />
                                   
                                </Routes>
                            </StyledSite>
                        </BrowserRouter>
                        <ErrorDialog />
                    </AccessTokenContext.Provider>
                </>
            )}
        </>
    );
}

