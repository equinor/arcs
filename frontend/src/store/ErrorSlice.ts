import { createSlice, PayloadAction } from "@reduxjs/toolkit";

interface ErrorState {
    errorMessage: string;
    isError: boolean;
}

const initialState: ErrorState = {
    errorMessage: "",
    isError: false,
};

export const errorSlice = createSlice({
    name: "error",
    initialState,
    reducers: {
        showError: (state, action: PayloadAction<string>) => {
            state.errorMessage = action.payload;
            state.isError = true;
        },
        dismissError: (state) => {
            state.errorMessage = "";
            state.isError = false;
        },
    },
});

export const { showError: showError, dismissError: dismissError } = errorSlice.actions;
export default errorSlice.reducer;
