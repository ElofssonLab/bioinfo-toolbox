Here are the commands needed to get GDT_TS values as in CASP:
=============================================================

fix_numbering.pl model.pdb native.pdb >model.pdb.fixed

```
cd <LGA_package_path>
ulimit -s unlimited
./runlga.mol_mol.pl full/path/to/model.pdb.fixed full/path/to/native.pdb -3  -ie  -o1  -sda  -d:4  -gdc_sc  -swap -ch1:A -ch2:A

gdt=`grep 'SUMMARY(GDT)' ./RESULTS/model.pdb.fixed.target.pdb.res | awk '{print $7}'`
```


Comments:
=========

You need to run fix_numbering.pl script prior to running LGA tool. You can find the script in Sscore directory.

runlga.mol_mol.pl is a wrapper script written by the package author himself. It puts the files in the correct format that is needed to run "lga" binary.

Options "-3  -ie  -o1  -sda  -d:4  -gdc_sc  -swap" are exactly the options that they use in CASP.

Obviously, you need to change -ch1:A -ch2:A parameters to match the chains in your proteins. If your protein does not have chain labeled, then leave out the corresponding parameter. This step is CRUCIAL, otherwise program crashes.

The results will always be in RESULTS directory of LGA_package. Unfortunately, there is no way to redirect the results to a different directory (or I just didn't find one). The output file is ./RESULTS/model.pdb.fixed.target.pdb.res. The grep line that I wrote extracts the GDT_TS score.
