# DoSM

[![CI](https://github.com/yui0303/Matrix-Diagonalization/actions/workflows/tests.yml/badge.svg)](https://github.com/yui0303/Matrix-Diagonalization/actions/workflows/tests.yml)

## Features

- Implement QR decomposition using two approaches, Gram-Schmidt and Householder.
- Pure c++
- Wrapping the implementation of c++ into python using Pybind11

## Build the Environment

- Clone the repo
- To install the these minimal pre-requisites on Ubuntu/Debian like linux operating systems, execute (in a terminal):
    ```
    sudo apt-get update
    sudo apt-get install -y build-essential make python3 python3-pip python3-pytest python3-numpy python3-scipy
    ```
- Install the Python dependencies
    ```
    python3 -m pip install --upgrade pip
    python3 -m pip install -r python/requirements.txt
    ```
- Build the system
    ```
    make
    ```

## Getting start

- Try `make test` to validate the correctness of the system
- Try `make demo` to run the benchmark
- Export to `PYTHONPATH`
    ```bash
    export MATRIX_PATH=/path/to/cpp
    export PYTHONPATH=${MATRIX_PATH}:${PYTHONPATH}
    ```

## User Tutorial
[API Introduction](https://hackmd.io/ohdBhIx7QGStMvYamwVtxg)