import { createSlice } from "@reduxjs/toolkit";
import type { PayloadAction } from "@reduxjs/toolkit";

type AuthState = {
    username: string;
    token: string;
    environment: string;
    homeAccountId: string;
    tenantId: string;
    localAccountId: string;
};

const authSlice = createSlice({
    name: "auth",
    initialState: {
        username: "",
        token: "",
        environment: "",
        homeAccountId: "",
        tenantId: "",
        localAccountId: "",
    } as AuthState,
    reducers: {
        setCredentials: (
            state,
            {
                payload: { username, token, environment, homeAccountId, tenantId, localAccountId },
            }: PayloadAction<{
                username: string;
                token: string;
                environment: string;
                homeAccountId: string;
                tenantId: string;
                localAccountId: string;
            }>
        ) => {
            state.username = username;
            state.token = token;
            state.environment = environment;
            state.homeAccountId = homeAccountId;
            state.tenantId = tenantId;
            state.localAccountId = localAccountId;
        },
    },
});

export const { setCredentials } = authSlice.actions;

export default authSlice.reducer;
