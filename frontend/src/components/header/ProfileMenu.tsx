import { useMsal } from "@azure/msal-react";
import { Button, Icon, Menu, Tooltip, Typography } from "@equinor/eds-core-react";
import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { Icons } from "../../utils/icons";

const ProfileMenu = () => {
    const [isOpen, setIsOpen] = useState<boolean>(false);
    const navigate = useNavigate();
    const { instance, accounts } = useMsal();
    const [anchorEl, setAnchorEl] = useState<HTMLButtonElement | null>(null);

    const handleLogOut = () => {
        instance.logout({
            onRedirectNavigate: () => {
                return false;
            },
        });
        navigate("/");
        closeMenu();
    };

    const openMenu = () => {
        setIsOpen(true);
    };

    const closeMenu = () => {
        setIsOpen(false);
    };

    const getName = () => {
        return accounts[0] === undefined ? "" : accounts[0].name;
    };

    const getUserName = () => {
        return accounts[0] === undefined ? "" : accounts[0].username;
    };

    return (
        <>
            <Tooltip title="Account">
                <Button
                    ref={setAnchorEl}
                    id="anchor-default"
                    aria-haspopup="true"
                    aria-expanded={isOpen}
                    aria-controls="menu-default"
                    onClick={() => (isOpen ? closeMenu() : openMenu())}
                >
                    <Icon name={Icons.Account} />
                </Button>
            </Tooltip>
            <Menu
                open={isOpen}
                id="menu-default"
                aria-labelledby="anchor-default"
                onClose={closeMenu}
                anchorEl={anchorEl}
                style={{ width: "275px" }}
            >
                <Menu.Item>
                    <Typography style={{ display: "block" }} variant="body_short" group="paragraph">
                        {getName()}
                    </Typography>
                    <Typography variant="overline" group="paragraph">
                        {getUserName()}
                    </Typography>
                </Menu.Item>
                <Menu.Item onClick={handleLogOut}>
                    <Icon name={Icons.LogOut} />
                    Sign out
                </Menu.Item>
            </Menu>
        </>
    );
};

export default ProfileMenu;
