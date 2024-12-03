import logging
import os
from typing import Annotated, Any

import httpx
import jwt
import jwt.algorithms
from dotenv import load_dotenv
from fastapi import HTTPException, Security
from fastapi.security import OAuth2AuthorizationCodeBearer

load_dotenv()

logger = logging.getLogger(__name__)

OPEN_ID_CONFIG_URI = "https://login.microsoftonline.com/3aa4a235-b6e2-48d5-9195-7fcf05b459b0/v2.0/.well-known/openid-configuration"
CLIENT_ID = os.environ.get("CLIENT_ID", "")


def _fetch_openid_configuration(auth_endpoint: str) -> Any:
    oid_conf_response = httpx.get(auth_endpoint)
    oid_conf_response.raise_for_status()
    return oid_conf_response.json()


oid_conf = _fetch_openid_configuration(OPEN_ID_CONFIG_URI)

jwks_client = jwt.PyJWKClient(oid_conf["jwks_uri"])

oauth2_scheme = Security(
    OAuth2AuthorizationCodeBearer(
        authorizationUrl=oid_conf["authorization_endpoint"],
        tokenUrl=oid_conf["token_endpoint"],
        auto_error=False,
        scopes={os.environ.get("API_SCOPE", ""): "Access to the arcs API"},
    )
)


def authenticated_user_claims(
    jwt_token: Annotated[str, oauth2_scheme],
) -> Any:
    if not jwt_token:
        raise HTTPException(401, "Missing token in Authorization header")
    try:
        signing_key = jwks_client.get_signing_key(
            jwt.get_unverified_header(jwt_token)["kid"]
        )
        claims = jwt.decode(
            jwt_token,
            key=signing_key,
            algorithms=["RS256"],
            audience=["api://" + CLIENT_ID, CLIENT_ID],
        )
        return claims
    except jwt.exceptions.InvalidTokenError as e:
        raise HTTPException(401, str(e))
