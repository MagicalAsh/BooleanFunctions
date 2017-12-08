coverage run testBooleanTools.py
coverage html
rm -r ~/public_html/documentation/htmlcov/
mv htmlcov/ ~/public_html/documentation/
chmod -R 755 ~/public_html/
