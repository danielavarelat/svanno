# svanno
Annotate SVs with breakpoint intersecting exons/introns/UTRs of OncoKB oncogenes/tumor suppressors

# setup
`pip install -r requirements.txt`

# usage
### comandline environment
`python svanno.py -i <vcf> -o <bedpe>`
### python environment
```
import sys
sys.path.append(<path to this repo>)
import svanno
import vcf2bedpe
bedpe = svanno.annotate(vcf2bedpe.convert(<path to vcf file>)) 
```

# example output
```
chr1	start1	end1	chr2	start2	end2	name	score	strand1	strand2	type	annot
1	2000	2001	1	3000	3001	DUP(1:2001-1:3001)	100	+	+	DUP	Intergenic>>SDHB:exon8
2	400 401	2	500	501	DEL(2:401-2:501)	100	+	+	DEL	ASXL2:intron2>>ASXL2:5pUTR
```

![Alt text](sv_annot_meme.png?raw=true "")
