
for NB in notebooks:
	jupyter nbconvert --to notebook --inplace --execute $NB
