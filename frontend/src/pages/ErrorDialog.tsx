import { Button, Dialog, Typography } from "@equinor/eds-core-react";
import { useAppDispatch, useAppSelector } from "../store/Hooks";
import { dismissError } from "../store/ErrorSlice";

const ErrorDialog = () => {
    const error = useAppSelector((state) => state.error);
    const dispatch = useAppDispatch();

    const onDismissError = () => {
        dispatch(dismissError());
    };

    return (
        <Dialog open={error.isError} style={{ width: "30rem" }}>
            <Dialog.Header>
                <Dialog.Title>Error</Dialog.Title>
            </Dialog.Header>
            <Dialog.CustomContent scrollable={true}>
                <Typography style={{ marginBottom: "15px" }}>{error.errorMessage}</Typography>
            </Dialog.CustomContent>
            <Dialog.Actions>
                <Button onClick={onDismissError}>Dismiss Error</Button>
            </Dialog.Actions>
        </Dialog>
    );
};

export default ErrorDialog;
