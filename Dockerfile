FROM python:3.12-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY pmidfrag ./pmidfrag
COPY tests ./tests
RUN pip install --no-cache-dir . pytest && python -m pytest -q tests
ENTRYPOINT ["pmidfrag"]
