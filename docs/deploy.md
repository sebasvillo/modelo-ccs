# Deploying the API

The API runs as a Docker container on Render's free plan (`render.yaml`, `Dockerfile`).
Free instances sleep after 15 min without traffic; the first request then takes about a minute.

- Every push to `main` rebuilds and redeploys (`autoDeploy: true`).
- CI builds the same image and smoke-tests `/health` and `/v1/case` (`.github/workflows/ci.yml`, job `docker`).
- CORS origins come from `CCS_CORS_ORIGINS` (comma-separated) in `render.yaml`.
- Custom domain: `api.sebasvillo.com` → CNAME to the service's `onrender.com` host, DNS-only
  (not proxied) so Render can issue the TLS certificate.
