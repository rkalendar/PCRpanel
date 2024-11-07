## PCRpanel 
## "Custom Amplicon Panels Designer is a tool for Amplicon Sequencing Technology, an ultra-high multiplex tiling PCR sequencing approach targeted DNA sequencing technology for next-generation sequencing (NGS) and nanopore sequencing (ONT)."
by Ruslan Kalendar 

email: ruslan.kalendar@helsinki.fi

[Web](https://primerdigital.com/tools/)

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 23 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

```
java -jar <PCRpanelPath>\dist\PCRpanel.jar <PCRpanelPath>\test\config.file

java -jar C:\PCRpanel\dist\PCRpanel.jar C:\PCRpanel\test\config.file 
```

### Basic usage: 
To enter parameters and specify the location of the target files and primer's file, you must specify this via a file on the command line. An example of such a file here (file name or extension does not matter):

> **config.file**
```
target_path=C:\MyPrograms\Java\PCRpanel\test\NG_012059.txt
target_primers=C:\MyPrograms\Java\PCRpanel\test\primers.txt
minPCR=100
maxPCR=500
minLen=18
maxLen=24
minTm=60
maxTm=63
minLC=80
3end=w
5end=
forwardtail=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
reversetail=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

```


## The output is saved in tab-delimited, plain text files. 




