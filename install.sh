echo "Installing eggnogmapper environment"
conda env create -f envs/eggnogmapper.yaml
echo "Installing metasequencing environment (this may take some time)"
conda env create -f envs/metasequencing.yaml

