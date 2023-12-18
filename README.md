# DoSM

[![CI / Ubuntu-22.04](https://github.com/yui0303/Matrix-Diagonalization/actions/workflows/tests.yml/badge.svg)](https://github.com/yui0303/Matrix-Diagonalization/actions/workflows/tests.yml)

## Build the Environment

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

- Build the project
```
make
```

## How to use

- Test the project
```
make test
```
- Demo the project
```
make demo
```
