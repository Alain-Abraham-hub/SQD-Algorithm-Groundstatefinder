# Docker Setup for SQD-Algorithm-Groundstatefinder

This guide explains how to use Docker to run the SQD Algorithm Groundstatefinder project.

## Prerequisites

- [Docker](https://www.docker.com/products/docker-desktop)
- [Docker Compose](https://docs.docker.com/compose/) (optional, but recommended)

## Quick Start (Pre-built Image)

The image is already built and available on Docker Hub, so you can skip building and use it directly:

```bash
# Pull the pre-built image
docker pull alain1435/sqd-groundstatefinder:latest
```

## Building the Docker Image Locally

If you want to build the image yourself:

```bash
# Build the image
docker build -t sqd-groundstatefinder .
```

## Running with Docker

All examples below use the pre-built image `alain1435/sqd-groundstatefinder:latest`. If you built a local image, replace this with your image name.

### Interactive Python Shell

```bash
docker run -it --rm -v $(pwd):/app alain1435/sqd-groundstatefinder:latest python
```

### Running Python Scripts

```bash
# Run a specific Python script
docker run --rm -v $(pwd):/app alain1435/sqd-groundstatefinder:latest python visualize_circuit.py

# Run with IBM Quantum credentials (if needed)
docker run --rm \
  -v $(pwd):/app \
  -e QISKIT_IBM_TOKEN="your_token" \
  -e QISKIT_IBM_INSTANCE="ibm-q/open/main" \
  alain1435/sqd-groundstatefinder:latest \
  python -c "from src.ansatz.run_on_ibm import run_ansatz_on_ibm; ..."
```

### Running Jupyter Notebooks

```bash
# Start Jupyter server
docker run --rm \
  -p 8888:8888 \
  -v $(pwd):/app \
  alain1435/sqd-groundstatefinder:latest \
  jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

Then access at `http://localhost:8888`

## Using Docker Compose

The `docker-compose.yml` file is configured to use the pre-built image from Docker Hub.

### Start Services

```bash
# Start the main service
docker-compose run --rm sqd-groundstatefinder

# Start Jupyter server
docker-compose up jupyter
```

### Run Commands

```bash
# Run a Python script
docker-compose run --rm sqd-groundstatefinder python visualize_circuit.py

# Open interactive shell
docker-compose run --rm sqd-groundstatefinder python
```

### Stop Services

```bash
docker-compose down
```

## Mount Volumes

All examples above mount the current directory (`$(pwd)`) to `/app` inside the container, allowing live editing of code and persistence of results.

## Environment Variables

For IBM Quantum integration, set environment variables:

```bash
docker run --rm \
  -v $(pwd):/app \
  -e QISKIT_IBM_TOKEN="your_token" \
  -e QISKIT_IBM_INSTANCE="ibm-q/open/main" \
  alain1435/sqd-groundstatefinder:latest \
  python your_script.py
```

Or define them in a `.env` file and use `--env-file`:

```bash
# Create .env file in project root
echo "QISKIT_IBM_TOKEN=your_token" >> .env
echo "QISKIT_IBM_INSTANCE=ibm-q/open/main" >> .env

docker run --rm \
  -v $(pwd):/app \
  --env-file .env \
  alain1435/sqd-groundstatefinder:latest \
  python your_script.py
```

## Troubleshooting

### Permission Issues with Jupyter

If you get permission errors with Jupyter, try:

```bash
docker run --rm \
  -p 8888:8888 \
  -v $(pwd):/app \
  -u root \
  alain1435/sqd-groundstatefinder:latest \
  jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

### Large Dependencies

The image includes all scientific computing dependencies (OpenBLAS, LAPACK, etc.) which may take time on first build. Subsequent builds will be faster due to Docker layer caching.

### GPU Support (Optional)

For GPU acceleration with Qiskit, you would need to build the image locally with a GPU-enabled base image:

```dockerfile
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04
# ... rest of Dockerfile
```

And run with:

```bash
docker run --rm --gpus all -v $(pwd):/app your-custom-image python
```
