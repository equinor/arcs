import { ThunkAction, UnknownAction, configureStore } from "@reduxjs/toolkit";
import progressReducer from "./ProgressSlice";
import { apiSlice } from "./Api/ApiSlice";
import errorHandlerMiddleware from "./ErrorHandlerMiddleware";
import authReducer from "./AuthSlice";
import errorReducer from "./ErrorSlice";

const store = configureStore({
    reducer: {
        progress: progressReducer,
        auth: authReducer,
        error: errorReducer,
        [apiSlice.reducerPath]: apiSlice.reducer,
    },
    middleware: (getDefaultMiddleware) => getDefaultMiddleware().concat(apiSlice.middleware, errorHandlerMiddleware),
});

export type RootState = ReturnType<typeof store.getState>;
export type AppDispatch = typeof store.dispatch;
export type AppThunk<ReturnType = void> = ThunkAction<ReturnType, RootState, unknown, UnknownAction>;
export default store;
