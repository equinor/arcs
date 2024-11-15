import { Button, Icon, Tooltip, TopBar, Typography } from "@equinor/eds-core-react";
import styled from "styled-components";
import { Icons } from "../../utils/icons";
import { useState } from "react";
import config from "../../configuration";
import ProfileMenu from "./ProfileMenu";
import { useAppSelector } from "../../store/Hooks";
import HelpDialog from "../../pages/Help/HelpDialog";



import { ReactSVG } from "react-svg";

const StyledCustomContent = styled.div`
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    align-items: flex-start;
`;

const StyledTitle = styled.div`
    display: flex;
    align-items: center;
    gap: 0.8rem;
`;

const StyledTitleText = styled.div`
    display: flex;
    align-items: left;
    flex-direction: column;
`;

const IconStyle = styled.div`
    display: flex;
    align-items: center;
    flex-direction: row-reverse;
    > * {
        margin-left: 1rem;
    }
`;

const HandPointer = styled.div`
    cursor: pointer;
`;

export function Header() {
    const [isHelpOpen, setIsHelpOpen] = useState(false);

    const handleHelpClick = () => {
        setIsHelpOpen(true);
    };

    const handleCloseHelp = () => {
        setIsHelpOpen(false);
    };

    return (
        <TopBar>
            <HandPointer>
                <TopBar.Header
                    onClick={() => {
                        window.location.href = "/";
                    }}
                >
                    <StyledTitle>
                        <StyledTitleText>
                             <Typography variant="body_long_bold" color="primary">
                                CO2 impurities tool
                            </Typography>
                        </StyledTitleText>
                    </StyledTitle>
                </TopBar.Header>
            </HandPointer>
            <TopBar.Actions>
                <IconStyle>
                    <Button variant="ghost" onClick={() => window.open(config.FEEDBACK_FORM_URL)}>
                        <Typography>Give Feedback</Typography>
                        <Icon name={Icons.ExternalLink} />
                    </Button>
                </IconStyle>
                <Tooltip title="Help">
                    <Button onClick={handleHelpClick} style={{ marginRight: "10px" }}>
                        <Icon name={Icons.HelpOutline} />
                    </Button>
                </Tooltip>
                <ProfileMenu />
            </TopBar.Actions>
            <HelpDialog isHelpOpen={isHelpOpen} handleCancel={handleCloseHelp} />
        </TopBar>
    );
}
