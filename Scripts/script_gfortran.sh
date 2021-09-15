#!/bin/bash
#echo `pwd`
#escolher diretorio
echo -n "Introduz diretório:"
read dir1



#default para o diretorio de um file, para nao ser tão chato
if [ $dir1 == "d" ]; then
	dir1="/home/$USER/Documents/Estágio/TOV/"
	echo "Default: $dir1"
fi


#ir para o diretorio escolhido
cd $dir1

#pwd

#nome d0 ficheiro
echo -n "Nome do ficheiro:"
read file

#versao d0 fotran, so numeros
echo -n "Versão de Fortran:"
read version

#echo "$version"

#escolher versao default
default=90

#default para a versao de fortran, para nao ser tão chato
if [ $version == "d" ]; then
	version="$default"
	echo "Default: f$default"
fi

#introduz um f na versão,
version="f$version"

#echo "$version"

#chama o compilador
echo `gfortran -o $file $file.$version`
echo "-------------- Programa OutPut ---------------"
#corre o script
echo `./$file`
