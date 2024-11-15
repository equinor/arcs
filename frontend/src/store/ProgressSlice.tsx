import { createSlice, PayloadAction } from "@reduxjs/toolkit";


interface ProgressState {
    runCalculation: Record<string, number>;
    autoGenerateCable: Record<string, number>;
}
enum MessageTypes {
    warning,
    error,
    info,
    progress_calculation,
    progress_autogenerate,
    progress_pdf,
}
interface Message {
    scenarioId: string;
    type: MessageTypes;
    value: string;
    datetime: string | undefined;
}

const initialState: ProgressState = {
    runCalculation: {} as Record<string, number>,
    autoGenerateCable: {} as Record<string, number>,
};

export const progressSlice = createSlice({
    name: "progress",
    initialState,
    reducers: {
        updateCalculationProgress: (state, action: PayloadAction<Message>) => {
            state.runCalculation[action.payload.scenarioId] = Number(action.payload.value);
            if (state.runCalculation[action.payload.scenarioId] > 100) {
                state.runCalculation[action.payload.scenarioId] = 100;
            }
        },
        
    },
});

export const {
    updateCalculationProgress: updateCalculationProgress,
} = progressSlice.actions;
export default progressSlice.reducer;
