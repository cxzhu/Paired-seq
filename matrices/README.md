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

For the Adult Brain cells only matrix, the “**Cluster**” column in “**Cell_embeddings_AdultBrain_Only.xls**” is the cell type group annotations as below:

|Cluster|Annotation|
|-----|-----|
|1|AS|
|2|MG|
|3|OC|
|4|Ex1|
|5|Ex2|
|6|Ex3|
|7|In1|
|8|In2|
|9|In3|



For the Fetal and Adult cells mixed matrix, the “**Ident**” column in “**Forebrain_Adult_Mix_Embeddings.xls**” is the cell type group annotations as below:

|Ident|Annotation|
|-----|-----|
|1|NP|
|2|dEx1|
|3|dEx2|
|4|dEx3|
|5|dIn1|
|6|dIn2|
|7|dIn3|
|8|dAS|
|9|Unknown|
