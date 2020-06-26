echo "Installing python2 dependent environment"
conda env create -f envs/metasequencing_py2.yaml
echo "Installing main metasequencing environment (this may take some time)"
conda env create -f envs/metasequencing.yaml

