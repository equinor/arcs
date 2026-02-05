FROM docker.io/library/python:3.12-slim
RUN apt-get update -y && apt-get install -y gcc
RUN python -m pip install --upgrade pip
RUN pip install --upgrade setuptools

WORKDIR /app
COPY . .
RUN python -m pip install .


EXPOSE 5001
USER 1000
ENTRYPOINT uvicorn api.app:app --host 0.0.0.0 --port 5001