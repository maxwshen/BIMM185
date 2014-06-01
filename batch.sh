for file in ~/bimm185/BIMM185/pam/*
do
  echo $file
  /usr/bin/python2.6 localAlign.py inaa.fasta $file
done