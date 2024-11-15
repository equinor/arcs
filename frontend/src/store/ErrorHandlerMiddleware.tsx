import { isRejectedWithValue } from "@reduxjs/toolkit";
import type { MiddlewareAPI, Middleware } from "@reduxjs/toolkit";
import { showError } from "./ErrorSlice";

const errorHandlerMiddleware: Middleware = (api: MiddlewareAPI) => (next) => (action) => {
    if (isRejectedWithValue(action)) {
        const payload = action.payload as any;
        if (payloadHasDescription(payload)) {
            if (isTokenExpired(payload)) return next(action);
            api.dispatch(showError(payload.data.description));
        }
    }

    return next(action);
};

const payloadHasDescription = (payload: any) => {
    return payload.data !== undefined && payload.data.description !== undefined;
};

const isTokenExpired = (payload: any) => {
    return payload.status === 401 && payload.data.description.includes("expired");
};

export default errorHandlerMiddleware;
