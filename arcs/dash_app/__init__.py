from .server import app


# Re-export internal Flask server for consumption by gunicorn
server = app.server


__all__ = [
    "app",
    "server",
]
