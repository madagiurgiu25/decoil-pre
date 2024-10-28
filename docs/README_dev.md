# Build and release Decoil

```
cd decoil-pre
rm -rf decoil.egg-info/ build/ dist/ sdist/
python -m pip install setuptools wheel twine 
python -m pip install --upgrade pip
python setup.py sdist bdist_wheel
python -m twine upload --repository testpypi --verbose dist/*


```
