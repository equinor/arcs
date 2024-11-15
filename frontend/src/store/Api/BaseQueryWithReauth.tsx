import { fetchBaseQuery } from "@reduxjs/toolkit/query";
import type { BaseQueryFn, FetchArgs, FetchBaseQueryError } from "@reduxjs/toolkit/query";
import { fetchAccessToken } from "../../api/AuthConfig";
import { setCredentials } from "../AuthSlice";
import config from "../../configuration";
import { msalInstance } from "../../main";

export const baseQueryWithReauth: BaseQueryFn<string | FetchArgs, unknown, FetchBaseQueryError> = async (
    args,
    api,
    extraOptions
) => {
    const state = api.getState() as any;
    const token = state.auth.token;
    const environment = state.auth.environment;
    const homeAccountId = state.auth.homeAccountId;
    const username = state.auth.username;
    const tenantId = state.auth.tenantId;
    const localAccountId = state.auth.localAccountId;

    const instance = msalInstance;
    const baseQuery = fetchBaseQuery({
        baseUrl: `${config.BACKEND_URL}`,
        prepareHeaders: (headers) => {
            if (token) {
                headers.set("authorization", `Bearer ${token}`);
            }
            headers.set("Content-Type", "application/json");
            return headers;
        },
    });

    let result = await baseQuery(args, api, extraOptions);
    if (result.error && result.error.status === 401) {
        const refreshedToken = await fetchAccessToken(
            instance,
            environment,
            homeAccountId,
            username,
            tenantId,
            localAccountId
        );
        api.dispatch(
            setCredentials({
                username: username,
                token: refreshedToken,
                environment: environment,
                homeAccountId: homeAccountId,
                tenantId: tenantId,
                localAccountId: localAccountId,
            })
        );
        const baseQueryWithRefreshedToken = fetchBaseQuery({
            baseUrl: `${config.BACKEND_URL}`,
            prepareHeaders: (headers) => {
                if (refreshedToken) {
                    headers.set("authorization", `Bearer ${refreshedToken}`);
                }
                headers.set("Content-Type", "application/json");
                return headers;
            },
        });
        result = await baseQueryWithRefreshedToken(args, api, extraOptions);
    }
    return result;
};
