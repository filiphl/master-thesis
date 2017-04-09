for f in Surface*.dump
do
mkdir ${f:0:${#f}-5}
mv $f ${f:0:${#f}-5}/
done
