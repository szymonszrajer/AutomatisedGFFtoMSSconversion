# Automatised GFF to MSS conversion
--final script and documentation are under construction--
file conversion inspired by GFF2MSS by Taro Maeda


Most, if not every annotation software used nowadays generates GFF files to describe genes and other features of sequences. This creates an issue, as none of the INSDC databases (DDBJ, ENA, NCBI) accepts annotations submitted in that format yet.

In the project that resulted in this program, we chose DDBJ for a submission of the *Gryllus Bimacaulatus* genome annotation.  [...]


To perform the first step of file preparation - format conversion, we tried using three scripts available online, all of which failed to facilitate the file conversion of an animal genome.

Apart from converting the file format, the errors that are a result of an automatic annotation also need to be resolved before the submission. To look for the errors within the annotation file, we used two validation programs distributed by the DDBJ - transChecker and jParser (https://www.ddbj.nig.ac.jp/ddbj/ume-e.html).

The Automatised GFF to MSS conversion program resolves some of the most common errors 
