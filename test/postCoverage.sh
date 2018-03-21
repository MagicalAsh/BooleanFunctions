coverage3 run testBooleanTools.py
coverage3 html
rm -r ~/public_html/documentation/htmlcov/
mv htmlcov/ ~/public_html/documentation/
chmod -R 755 ~/public_html/
