# API image: Python 3.13 + locked dependencies, no dev tools.
FROM python:3.13-slim

COPY --from=ghcr.io/astral-sh/uv:0.12 /uv /bin/uv
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy UV_PYTHON_DOWNLOADS=never PYTHONUNBUFFERED=1

WORKDIR /app
COPY pyproject.toml uv.lock README.md LICENSE ./
COPY src ./src
COPY api ./api
RUN uv sync --locked --no-dev

RUN useradd --create-home app
USER app
ENV PATH="/app/.venv/bin:$PATH"

# Render (and most hosts) pass the port in $PORT.
CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8000} --proxy-headers"]
