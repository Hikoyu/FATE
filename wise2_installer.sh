#!/bin/bash

# clangを確認
if [ ! `which clang` ] ; then
	echo "clang not installed"
	exit 1
fi

# pkg-configを確認
if [ ! `which pkg-config` ] ; then
	echo "pkg-config not installed"
	exit 1
fi

# glib2を確認
pkg-config glib-2.0
if [ $? == 1 ] ; then
	echo "glib2 not installed"
	exit 1
fi

# ファイル書き換え
sed 's/getline/getline_new/' src/HMMer2/sqio.c  > src/HMMer2/sqio.c~
mv src/HMMer2/sqio.c~ src/HMMer2/sqio.c

sed 's/isnumber/isdigit/' src/models/phasemodel.c > src/models/phasemodel.c~
mv src/models/phasemodel.c~ src/models/phasemodel.c

sed 's/csh welcome.csh/sh welcome.csh/' src/makefile > src/makefile~
mv src/makefile~ src/makefile

file_list=`find . -name "makefile"`
for file in $file_list; do
	sed -e 's/CC = cc/CC = clang/g' -e 's/glib-config/pkg-config glib-2.0/g' $file > "$file~"
	mv "$file~" $file
done

# コンパイル
cd src && make all
