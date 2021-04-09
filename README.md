# Landscape and Selection of Vaccine Epitopes in SARS-CoV-2

Preprint submitted to biorxiv: https://www.biorxiv.org/content/10.1101/2020.06.04.135004v1

Several larger files which were used part of these analyses are included alongside supplemental data on Mendeley due to GitHub file size limits: https://data.mendeley.com/datasets/c6pdfrwxgj/3.

## Instructions

All code used to analyze data and generate figures from this manuscript have been released within this github page, with some larger files (>100mb) stored on Mendeley data (https://data.mendeley.com/datasets/c6pdfrwxgj/3). In order to replicate analysis presented in this manuscript, follow the steps listed below.  R v3.5.2 and Python v3.7.5 were used for these analyses:

1. Clone this repository onto your local computer
2. Download the large data files (along with supplemental tables) from the above Mendeley data link.
3. Unpack the data files from above.  The subdirectories contained in the "Large_files" directory correspond to the same order as those within the github repository.  Merge these subdirectories together with your local repository files.
4. Run the following command in the repository root directory: 
   >tar -xvf ./SARS-CoV-2_epitope_landscape/Working/large_files_040521.tar.gz
5. For Figures 2-3, Fig. S1-S6, and Table S1-S8: Open the R file "Main_figures_resubmission.R" in the repository root directory, changing the "WORKING_ROOT" variable (line 2) to the path of the repository.  This file contains a step-by-step workflow for recreating these above listed figures.
6. Run `pip install -r requirements.txt`
7. Run `mhcflurry-downloads fetch models_class1_presentation`
8. To generate Figures 4-6, run `generate-python-figures.sh`, which runs a series of IPython notebooks to select B-cell epitope regions, T-cell epitopes, vaccine peptides, and the generation of all related plots and figures. 
 
## Python Package Versions

These package versions are also listed in `requirements.txt`:

```
absl-py==0.8.1
alabaster==0.7.12
appdirs==1.4.3
appnope==0.1.0
astor==0.7.1
astroid==2.3.1
astropy==3.0.4
astunparse==1.6.3
attrs==19.2.0
Babel==2.8.0
backcall==0.1.0
beautifulsoup4==4.9.1
biopython==1.72
bitarray==1.2.1
bleach==3.1.0
blis==0.4.1
Bottleneck==1.3.2
bsuite==0.3.2
bump2version==1.0.1
bumpversion==0.6.0
bz2file==0.98
CacheControl==0.12.6
cached-property==1.5.1
cachetools==4.0.0
catalogue==1.0.0
certifi==2019.11.28
chardet==3.0.4
chex==0.0.2
Click==7.0
cloudpickle==1.2.2
cmudict==0.4.3
commonmark==0.9.1
cycler==0.10.0
cymem==2.0.3
Cython==0.29.14
-e git+git@github.com:hammerlab/datacache.git@73bcac02d37cf153710a07fbdc636aa55cb214ca#egg=datacache
decorator==4.4.1
defusedxml==0.6.0
descartes==1.1.0
dill==0.3.2
dm-env==1.2
dm-tree==0.1.5
docopt==0.6.2
docutils==0.15.2
entrypoints==0.3
fastai==1.0.61
fastprogress==1.0.0
filelock==3.0.12
Flask==1.1.1
frozendict==1.2
funcsigs==1.0.2
future==0.18.2
gast==0.2.2
gevent==1.5a3
ghostscript==0.6
gin-config==0.3.0
google-api-python-client==1.7.11
google-auth==1.11.0
google-auth-httplib2==0.0.3
google-auth-oauthlib==0.4.1
google-pasta==0.2.0
googleapis-common-protos==1.52.0
greenlet==0.4.15
grpcio==1.31.0
gtfparse==1.1.2
gunicorn==20.0.4
gym==0.15.4
h5py==2.10.0
hdmedians==0.13
httplib2==0.17.0
idna==2.8
imageio==2.9.0
imagesize==1.2.0
importlib-metadata==0.23
importlib-resources==3.0.0
ipdb==0.12.3
ipykernel==5.1.2
ipython==7.11.1
ipython-genutils==0.2.0
ipywidgets==7.5.1
isort==4.3.21
itsdangerous==1.1.0
jax==0.1.75
jax-loss==0.0.1
jaxlib==0.1.52
jedi==0.16.0
Jinja2==2.11.2
joblib==0.14.1
jsonschema==3.1.1
jupyter==1.0.0
jupyter-client==5.3.4
jupyter-console==6.0.0
jupyter-core==4.6.0
Keras==2.3.1
Keras-Applications==1.0.8
Keras-Preprocessing==1.1.2
kiwisolver==1.1.0
lazy-object-proxy==1.4.2
lockfile==0.12.2
logomaker==0.8
lxml==4.5.0
m2r==0.2.1
Markdown==2.6.11
MarkupSafe==1.0
matplotlib==3.1.2
matplotlib-venn==0.11.5
mccabe==0.6.1
memoized-property==1.0.3
mesh-tensorflow==0.1.16
-e git+git@github.com:hammerlab/mhcflurry.git@3ee69c6d134338359542154a1107b0c447e1960a#egg=mhcflurry
-e git+git@github.com:til-unc/mhcgnomes.git@3e65405aefe7019f0bc751260ad73dd888e434ea#egg=mhcgnomes
-e git+git@github.com:iskandr/mhcnames.git@71694b9d620db68ceee44da1b8422ff436f15bd3#egg=mhcnames
-e git+git@github.com:hammerlab/mhctools.git@c64f5d9d08f2d474fe7868f42c8f116e532bcc35#egg=mhctools
mistune==0.8.4
mizani==0.7.1
mock==2.0.0
more-itertools==7.2.0
mpmath==1.1.0
msgpack==0.6.2
murmurhash==1.0.2
natsort==7.0.1
nbconvert==5.6.0
nbformat==4.4.0
networkx==2.4
ngs-plumbing==0.13.1
nltk==3.4.5
nose==1.3.7
notebook==6.0.1
np-utils==0.5.11.1
numexpr==2.7.1
numpy==1.18.1
nvidia-ml-py3==7.352.0
oauth2client==4.1.3
oauthlib==3.1.0
opencv-python==4.1.2.30
opt-einsum==3.1.0
packaging==20.1
palettable==3.3.0
pandas==1.0.0
pandocfilters==1.4.2
parso==0.6.0
patsy==0.5.1
pbr==4.2.0
pdfkit==0.6.1
pexpect==4.8.0
pickleshare==0.7.5
Pillow==7.2.0
pkginfo==1.5.0.1
plac==1.1.3
plotnine==0.7.1
pluggy==0.13.1
portalocker==2.0.0
preshed==3.0.2
progressbar33==2.4
prometheus-client==0.7.1
promise==2.3
prompt-toolkit==3.0.3
pronouncing==0.2.0
protobuf==3.13.0
psutil==5.6.3
ptpython==2.0.6
ptyprocess==0.6.0
py==1.8.1
pyasn1==0.4.8
pyasn1-modules==0.2.8
pycairo==1.18.2
pyfastx==0.6.5
pyflakes==2.0.0
pyglet==1.3.2
Pygments==2.5.2
pylint==2.4.2
pynvx==1.0.0
pyopa==0.8.0
pypandoc==1.4
pyparsing==2.4.6
Pyphen==0.9.5
pyrsistent==0.15.4
pysam==0.15.3
pytest==5.3.5
python-dateutil==2.8.1
python-Levenshtein==0.12.0
pytorch-lightning==0.9.0
pytypes==1.0b5
pytz==2019.3
PyVCF==0.6.8
PyWavelets==1.1.1
PyYAML==5.3.1
pyzmq==18.1.0
qtconsole==4.5.5
readme-renderer==24.0
recommonmark==0.6.0
regex==2020.7.14
requests==2.22.0
requests-oauthlib==1.3.0
requests-toolbelt==0.9.1
rlax==0.0.2
roman==3.0
rouge-score==0.0.4
rsa==4.0
sacrebleu==1.4.13
sacremoses==0.0.43
scikit-bio==0.5.5
scikit-image==0.17.2
scikit-learn==0.22.1
scipy==1.4.1
seaborn==0.11.1
Send2Trash==1.5.0
sentencepiece==0.1.91
seqlogo==5.29.7
sercol==0.1.4
setuptools-scm==3.5.0
shellinford==0.4.1
simplejson==3.16.0
six==1.14.0
snowballstemmer==2.0.0
soupsieve==2.0.1
spacy==2.3.2
Sphinx==2.4.0
sphinxcontrib-applehelp==1.0.1
sphinxcontrib-devhelp==1.0.1
sphinxcontrib-htmlhelp==1.0.2
sphinxcontrib-jsmath==1.0.1
sphinxcontrib-qthelp==1.0.2
sphinxcontrib-serializinghtml==1.1.3
srsly==1.0.2
statsmodels==0.12.0
sympy==1.5.1
t5==0.6.4
tensor2tensor==1.7.0
tensorboard==2.2.0
tensorboard-plugin-wit==1.7.0
tensorflow==2.3.0
tensorflow-datasets==3.2.1
tensorflow-estimator==2.3.0
tensorflow-metadata==0.23.0
tensorflow-text==2.3.0
termcolor==1.1.0
terminado==0.8.2
testpath==0.4.2
tfds-nightly==3.2.1.dev202008240105
thinc==7.4.1
tifffile==2020.8.25
tinytimer==0.0.0
tokenizers==0.8.1rc1
toolz==0.10.0
torch==1.4.0
torchvision==0.7.0
tornado==6.0.3
tqdm==4.48.2
traitlets==4.3.3
transformers==3.0.2
trax==1.3.4
twine==2.0.0
typechecks==0.1.0
typed-ast==1.4.0
Unidecode==1.1.2
uritemplate==3.0.1
urllib3==1.25.8
wasabi==0.7.1
wcwidth==0.1.8
webencodings==0.5.1
weblogo==3.7.1
Werkzeug==1.0.1
widgetsnbextension==3.5.1
wrapt==1.11.2
xlrd==1.1.0
XlsxWriter==1.1.0
xlwt==1.3.0
xvfbwrapper==0.2.9
zipp==0.6.0
```

