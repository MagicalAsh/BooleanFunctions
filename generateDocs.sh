#!/bin/bash
if [ ! -d ./docs ]; then
	mkdir docs
fi

pydocmd simple booleantools.BooleanFunction+ > docs/boolFunc.md
pydocmd simple booleantools+ > docs/boolTools.md

cd docs/

pydocmd new
echo "- BooleanFunction: boolean_functions.md << boolFunc.md" >> pydocmd.yml
echo "- BooleanTools: boolean_tools.md << boolTools.md" >> pydocmd.yml
sed -i "s/site_name.*/site_name: BooleanTools Documentation/g" pydocmd.yml

pydocmd build


