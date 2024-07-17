Compressed files are splitted due to filesize limite of 25MB.

```
### MacOS commandline
zip -s 25m -r Adult_Brain_split.zip ./Adult_Brain
zip -s 25m -r Fetal_adult_merged_split.zip ./Fetal_forebain_Adult_Brain_mix
```

To decompress the matrix files, download all splitted files into the same folder.
```
### MacOS commandline
zip -s 0 Adult_Brain_split.zip --out Adult_Brain_merged.zip
unzip Adult_Brain_merged.zip
zip -s 0 Fetal_adult_merged_split.zip --out Fetal_adult_merged_merged.zip
unzip Fetal_adult_merged_merged.zip
```
