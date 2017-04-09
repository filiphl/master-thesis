for f in Surface*.dump
do
mv ${f:0:${#f}-5}/* .
rm -rf ${f:0:${#f}-5}/
done
