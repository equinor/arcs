import { Button, Dialog } from "@equinor/eds-core-react";
import { FunctionComponent } from "react";
import HelpPage from "./HelpPage";

type Props = {
    isHelpOpen: boolean;
    handleCancel: Function;
};

const HelpDialog: FunctionComponent<Props> = ({ isHelpOpen, handleCancel }) => {
    const onClose = () => {
        handleCancel();
    };

    return (
        <Dialog open={isHelpOpen} isDismissable={true} onClose={onClose} style={{ width: "1000px", height: "85vh" }}>
            <Dialog.Title></Dialog.Title>
            <Dialog.CustomContent>
                <HelpPage />
            </Dialog.CustomContent>
            <Dialog.Actions>
                <Button onClick={onClose}>Close</Button>
            </Dialog.Actions>
        </Dialog>
    );
};

export default HelpDialog;
