# SV_caller
A python tool for identify lineage specific structure variation from multiple genome aligmnent
## Usage
first we extract maf into list format like this 
```
ref_chrom  ref_position sp1 sp2 sp3 sp4  ref2_pos
chr1  1  A  A  A  A  ref2_pos
chr1  2  T  T  A  T  ref2_pos
```
```
python 01.changemaf2list.py --maffile align.maf --speciesfile speciesfile --outfile out --sp1addlc ref2
```
here speciefile is file with each specie name per line which you want to extract.
ref2 is the reference speice name in outgroup to show the deletion in ingroup

then running

```
python 00.change2list.py --maffile maf.list --outgroup outgrouplist --out out.bed
```
here maf.list is the output list file from last step, outgroup list is the specie name per line for outgroups contrast to your target lineage.
output file is a bed like format, one SV per line, below in what each column represent
```
CHROM chromosome id of ref
START start location of SV
END end location of SV
START_O  start location for outgroup ref
END_O  end location for outgroup ref
IDENT_IN  sequence identity of SV for ingroup
IDENT_IN_UP  sequence identity upstream 50bp of SV for ingroup
IDENT_IN_DOWN  sequence identity downstream 50bp of SV for ingroup
IDENT_OUT_UP  sequence identity upstream 50bp of SV for outgroup
IDENT_OUT_DOWN  sequence identity upstream 50bp of SV for outgroup
LENGTH  length of SV
```
simply you can use awk to extract the SV meet your standard
```
awk '{if($6>=0.9 && $11 >= 5) print $0}' raw_sv.bed > filtered_sv.bed
```
you will get SVs over 5 bp and 90% sequences are identical for ingroup.
