# Build and release Decoil

```
cd decoil-pre

# install dependencies
python -m pip install setuptools wheel twine 
python -m pip install --upgrade pip

# build
rm -rf decoil.egg-info/ build/ dist/ sdist/
python -m pip install -r requirements.txt
python setup.py sdist bdist_wheel

# test installation locally from wheel
python -m pip install dist/decoil-*-py3-none-any.whl

# upload to testpypi
python -m twine upload --repository testpypi --verbose dist/*

# download from testpypi
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  decoil==<version>

```

# Debug

1. Error installing `parasails` on macos: 

````
ERROR: Could not build wheels for parasail, which is required to install pyproject.toml-based projects
```

That error usually happens on macOS when Python packages with C extensions (like parasail) fail to compile. parasail is a C library, and pip is trying to build it from source. Here’s how to fix it step by step:

```
xcode-select --install

# parasail requires some C libraries. On macOS, you can install dependencies via Homebrew
brew install gcc
brew install automake autoconf libtool

python3.10 -m pip install --upgrade pip setuptools wheel

python3.10 -m pip install parasails
```

2. Run individual test:

```
python -m pip install pytest
# -s will capture the print() and goes to console
pytest -s -k test_none_DR tests/encode/test_input.py 
```

3. Check coverage of the tests:

Overall summary:

```
pytest --cov=decoil
```

Check individual branch:

```
pytest --cov=decoil.encode.encode --cov-report=html --cov-branch
```